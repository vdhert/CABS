"""
Module for running cabsDock jobs.
"""

import re
import operator
from os import getcwd, mkdir
from os.path import exists, isdir, join, abspath
from time import sleep

from cabsDock.cluster import Clustering
from protein import ProteinComplex
from restraints import Restraints
from cabs import CabsRun
from utils import ProgressBar
from utils import SCModeler
from utils import plot_E_rmsds
from utils import plot_rmsd_N
from utils import check_peptide_sequence
from utils import _chunk_lst
from utils import mk_histos_series
from trajectory import Trajectory
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap
from filter import Filter

__all__ = ['Job']


class Job:
    """
    Class representing single cabsDock job.
    """

    def __repr__(self):
        return '\n'.join([k + ' : ' + str(v) for k, v in sorted(self.config.items())])

    def __init__(self, **kwargs):
        # TODO: defaults jesli job import
        self.config = kwargs

        # making sure work_dir is abspath
        self.config['work_dir'] = abspath(self.config['work_dir'])

        # checks if work_dir exists, creates it otherwise
        work_dir = self.config['work_dir']
        if exists(work_dir):
            if not isdir(work_dir):
                raise Exception('File %s already exists and is not a directory' % work_dir)
            # ans = raw_input('You are about to overwrite results in %s\nContinue? y or n: ' % work_dir)
            # if ans != 'y':
            #     exit(code=1)
        else:
            mkdir(work_dir)

        self.initial_complex = None

    def prepare_restraints(self):

        # generate receptor restraints
        receptor_restraints = Restraints(
            self.initial_complex.receptor.generate_restraints(*self.config['receptor_restraints'])
        )

        # additional restraints
        add_restraints = Restraints()

        if self.config['ca_rest_add']:
            add_restraints += Restraints(self.config['ca_rest_add'])

        if self.config['sg_rest_add']:
            add_restraints += Restraints(self.config['sg_rest_add'], sg=True)

        if self.config['ca_rest_file']:
            for filename in self.config['ca_rest_file']:
                add_restraints += Restraints(filename)

        if self.config['sg_rest_file']:
            for filename in self.config['sg_rest_file']:
                add_restraints += Restraints(filename, sg=True)

        receptor_restraints += add_restraints.update_id(self.initial_complex.new_ids)

        return receptor_restraints

    def process_excluding(self):
        pass

    def run_job(self):
        work_dir = self.config['work_dir']
        print('CABS-docking job {0}'.format(self.config['receptor']))

        # prepare initial complex
        print(' Building complex...')
        self.initial_complex = ProteinComplex(self.config)
        print(' ... done.')

        # run cabs
        print('CABS simulation starts.')
        cabs_run = CabsRun(self.initial_complex, self.prepare_restraints(), self.config)
        cabs_run.run()

        print('CABS simuation is DONE.')
        trajectory = cabs_run.get_trajectory()
        trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
        if self.config['native_pdb']:
            print('Calculating RMSD to the native structure...')
            print(
                'The native complex loaded from {0} consists of receptor (chain(s) {1}) and peptide(s) (chains(s) {2}).'
                .format(
                    self.config['native_pdb'],
                    self.config['native_receptor_chain'],
                    self.config['native_peptide_chain']
                )
            )
            rmslst = trajectory.rmsd_to_native(native_pdb=self.config['native_pdb'],
                                      native_receptor_chain=self.config['native_receptor_chain'],
                                      native_peptide_chain=self.config['native_peptide_chain'],
                                      model_peptide_chain=self.initial_complex.ligand_chains[0])
        elif self.config['reference_pdb']:
            rmslst = trajectory.rmsd_to_reference(
                                        ref_pdb=self.config['reference_pdb'],
                                        pept_chain=self.initial_complex.ligand_chains[0],
                                        align_mth=self.config['align']
                                        )
        trajectory.align_to(self.initial_complex.receptor)

        # energy fix
        trajectory.number_of_peptides = len(self.initial_complex.ligand_chains)
        tra, flt_inds = Filter(trajectory, N=1000).cabs_filter()
        tra.number_of_peptides = len(self.initial_complex.ligand_chains)

        rmsf_vals = _chunk_lst(trajectory.rmsf(self.initial_complex.receptor_chains), 15, 0)
        lbls = _chunk_lst([i.chid + str(i.resnum) + i.icode for i in trajectory.template.atoms if i.chid in self.initial_complex.receptor_chains], 15, "")

        pltdir = self.config['work_dir'] + '/plots'
        try: mkdir(pltdir)
        except OSError: pass

        mk_histos_series(rmsf_vals, lbls, pltdir + '/RMSF_target')

        if self.config['native_pdb'] or self.config['reference_pdb']:
            plot_E_rmsds(   [trajectory, tra],
                            [rmslst, rmslst[flt_inds,]],
                            ['total','interaction'],
                            pltdir + '/Ermsd')
            plot_rmsd_N(    rmslst.reshape(self.config['replicas'], -1),
                            pltdir + '/RMSDn')

        medoids, clusters_dict, clusters = Clustering(
            tra, 'chain ' + ','.join(self.initial_complex.ligand_chains)
        ).cabs_clustering()

        if self.config['contact_maps']:
            self.mk_cmaps(trajectory, medoids, clusters_dict, flt_inds, 4.5, pltdir)

        # Saving the trajectory to PDBs:
        trajectory.to_pdb(mode = 'replicas', to_dir = work_dir)
        # Saving top1000 models to PDB:
        tra.to_pdb(mode = 'replicas', to_dir= work_dir, name='top1000' )

        # Saving clusters in CA representation
        for i, cluster in enumerate(clusters):
            cluster.to_pdb(mode='replicas', to_dir=work_dir, name='cluster_{0}'.format(i))

        # Saving top10 models:
        if self.config['AA_rebuild']:
            # Saving top 10 models in AA representation:
            pdb_medoids = medoids.to_pdb()
            from cabsDock.ca2all import ca2all
            for i, file in enumerate(pdb_medoids):
                ca2all(file, output='model_{0}.pdb'.format(i), iterations=1, verbose=False)
        else:
            # Saving top 10 models in CA representation:
            medoids.to_pdb(mode='models', to_dir=work_dir, name='model')

        # dictionary holding results
        rmsds = [header.rmsd for header in medoids.headers]
        results = {}
        results['rmsds_10k'] = [header.rmsd for header in trajectory.headers]
        results['rmsds_1k'] = [header.rmsd for header in tra.headers]
        results['rmsds_10'] = rmsds
        results['lowest_10k'] = sorted(results['rmsds_10k'])[0]
        results['lowest_1k'] = sorted(results['rmsds_1k'])[0]
        results['lowest_10'] = sorted(results['rmsds_10'])[0]

        # Saving rmsd results
        with open(join(work_dir, 'rmsds.txt'), 'w') as outfile:
            outfile.write(str(results))

        # Not returned by default unless self.config['benchmark'] == True.
        if self.config['benchmark']:
            return results

    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, plots_dir):
        scmodeler = SCModeler(ca_traj.template)
        sc_traj_full = scmodeler.calculate_sc_traj(ca_traj.coordinates)
        sc_traj_1k = sc_traj_full.reshape(1, -1, len(ca_traj.template), 3)[:,top1k_inds,:,:]
        sc_med = scmodeler.calculate_sc_traj(meds.coordinates)
        shp = sc_med.shape
        sc_med = sc_med.reshape((shp[1], shp[0]) + shp[2:])

        cmapdir = self.config['work_dir'] + '/contact_maps'
        try: mkdir(cmapdir)
        except OSError: pass
        rchs = self.initial_complex.receptor_chains
        lchs = self.initial_complex.ligand_chains

        targ_cmf = ContactMapFactory(rchs, rchs, ca_traj.template)
        cmfs = {lig: ContactMapFactory(rchs, lig, ca_traj.template) for lig in lchs}
        cmap10ktarg = reduce(operator.add, targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg.zero_diagonal()
        cmap10ktarg.save_all(cmapdir + '/target_all')

        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(sc_traj_full, thr)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(cmapdir + '/replica_%i_ch_%s' % (n + 1, lig))
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + '/all_ch_%s' % lig)
            cmap10k.save_histo(plots_dir + '/all_contacts_histo_%s' % lig)
            cmap1k = cmf.mk_cmap(sc_traj_1k, thr)[0]
            cmap1k.save_all(cmapdir + '/top1000_ch_%s' % lig)
            cmaps_top = cmf.mk_cmap(sc_med, thr)
            for n, cmap in enumerate(cmaps_top):
                cmap.save_all(cmapdir + '/top_%i_ch_%s' % (n + 1, lig))
            for cn, clust in clusts.items():
                ccmap = cmf.mk_cmap(sc_traj_1k, thr, frames=clust)[0]
                ccmap.save_all(cmapdir + '/cluster_%i_ch_%s' % (cn, lig))

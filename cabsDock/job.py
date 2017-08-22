"""
Module for running cabsDock jobs.
"""

import operator
from os import getcwd, mkdir
from os.path import exists, isdir, abspath

from cabs import CabsRun
from cabsDock.cluster import Clustering
from cabsDock.cmap import ContactMapFactory
from cabsDock.plots import graph_RMSF
from cabsDock.plots import plot_E_RMSD
from cabsDock.plots import plot_RMSD_N
from cabsDock.utils import SCModeler
from filter import Filter
from protein import ProteinComplex
from restraints import Restraints
from trajectory import Trajectory

__all__ = ['Job']


class Job:
    """
    Class representing single cabsDock job.
    """

    def __repr__(self):
        return '\n'.join([k + ' : ' + str(v) for k, v in sorted(self.config.items())])

    def __init__(
            # TODO KEEP THE DEFAULTS UPDATED FOR TESTING.
            self,
            receptor=None,
            ligand=None,
            work_dir=getcwd(),
            replicas=10,
            mc_annealing=20,
            mc_cycles=50,
            mc_steps=50,
            temperature=(2.0, 1.0),
            replicas_dtemp=0.5,
            separation=20.0,
            insertion_clash=0.5,
            insertion_attempts=1000,
            ca_rest_add=None,
            ca_rest_file=None,
            ca_rest_weight=1.0,
            sc_rest_add=None,
            sc_rest_file=None,
            sc_rest_weight=1.0,
            receptor_restraints=('all', 4, 5.0, 15.0),
            dssp_command='mkdssp',
            fortran_command='gfortran',
            filtering_number=1000,  # number of models to filter
            filtering_fromeach=True,
            clustering_medoids=10,
            clustering_iterations=100,
            benchmark=False,
            AA_rebuild=True,
            contact_maps=True,
            reference_pdb=None,
            save_replicas=True,
            save_topn=True,
            save_clusters=True,
            save_medoids='AA',  # AA or CG. AA option requires MODELLER to be installed.
            load_cabs_files=None,
            align='SW',
            reference_alignment=None,
            save_config_file=True,
            image_file_format='svg',
            verbose=None,
            stride_command='stride',
            receptor_flexibility=None,
            exclude=None,
            save_cabs_files=None,
            output_clusters=False,
            peptide=None,
            add_peptide=None,
            random_seed=None,
            output_trajectories=None,
            config=None,
            no_aa_rebuild=False,
            excluding_distance=5.0,
            modeller_iterations=3,
            output_models=10
    ):
        if load_cabs_files and len(load_cabs_files) is 2:
            file_TRAF, file_SEQ = load_cabs_files
        else:
            file_TRAF = file_SEQ = None

        # TODO replace self.config dictionary with regular attributes. Clean up all usages of self.config in job methods.
        self.config = {
            'receptor': receptor,
            'ligand': ligand,
            'work_dir': work_dir,
            'replicas': replicas,
            'mc_cycles': mc_cycles,
            'mc_steps': mc_steps,
            'mc_annealing': mc_annealing,
            't_init': temperature[0],
            't_final': temperature[1],
            'replicas_dtemp': replicas_dtemp,
            'initial_separation': separation,
            'ligand_insertion_clash': insertion_clash,
            'ligand_insertion_attempts': insertion_attempts,
            'ca_rest_add': ca_rest_add,
            'ca_rest_file': ca_rest_file,
            'ca_rest_weight': ca_rest_weight,
            'sc_rest_add': sc_rest_add,
            'sc_rest_file': sc_rest_file,
            'sc_rest_weight': sc_rest_weight,
            'receptor_restraints': receptor_restraints,  # sequence gap, min length, max length
            'dssp_command': dssp_command,
            'fortran_compiler': fortran_command, # build (command, flags)
            'filtering': filtering_number,  # number of models to filter
            'filtering_fromeach': filtering_fromeach,
            'clustering_nmedoids': clustering_medoids,
            'clustering_niterations': clustering_iterations,  # number of clusters, iterations
            'benchmark': benchmark,
            'AA_rebuild': AA_rebuild,
            'contact_maps': contact_maps,
            'reference_pdb': reference_pdb,
            'save_replicas': save_replicas,
            'save_topn': save_topn,
            'save_clusters': save_clusters,
            'save_medoids': save_medoids,  # 'AA' or 'CG'. 'AA' option requires MODELLER to be installed.
            'file_TRAF': file_TRAF,
            'file_SEQ': file_SEQ,
            'align': align,
            'reference_alignment': reference_alignment,
            'save_config_file': save_config_file,
            'image_file_format': image_file_format,
            'receptor_flexibility': receptor_flexibility
        }

        # Job attributes collected.
        self.initial_complex = None
        self.restraints = None
        self.cabsrun = None
        self.trajectory = None
        self.filtered_trajectory = None
        self.filtered_ndx = None
        self.medoids = None
        self.clusters_dict = None
        self.clusters = None
        self.rmslst = {}
        self.results = None

        # Workdir processing:
        # making sure work_dir is abspath
        self.config['work_dir'] = abspath(self.config['work_dir'])

        # checks if work_dir exists, creates it otherwise
        work_dir = self.config['work_dir']
        if exists(work_dir):
            if not isdir(work_dir):
                raise Exception('File %s already exists and is not a directory' % work_dir)
        else:
            mkdir(work_dir)

    def prepare_restraints(self):

        # generate receptor restraints
        receptor_restraints = Restraints(
            self.initial_complex.receptor.generate_restraints(*self.config['receptor_restraints'])
        )

        # additional restraints
        add_restraints = Restraints('')

        if self.config['ca_rest_add']:
            add_restraints += Restraints.from_parser(self.config['ca_rest_add'])

        if self.config['sc_rest_add']:
            add_restraints += Restraints.from_parser(self.config['sc_rest_add'], sg=True)

        if self.config['ca_rest_file']:
            for filename in self.config['ca_rest_file']:
                add_restraints += Restraints.from_file(filename)

        if self.config['sc_rest_file']:
            for filename in self.config['sc_rest_file']:
                add_restraints += Restraints.from_file(filename, sg=True)

        receptor_restraints += add_restraints.update_id(self.initial_complex.new_ids)
        return receptor_restraints

    def cabsdock(self):
        ftraf = self.config.get('file_TRAF')
        fseq = self.config.get('file_SEQ')
        self.setup_job()
        withcabs = True if (ftraf is None or fseq is None) else False
        if withcabs:
            self.setup_cabs_run()
            self.execute_cabs_run()
        self.load_output(ftraf, fseq)
        self.score_results(n_filtered=self.config['filtering'], number_of_medoids=self.config['clustering_nmedoids'],
                           number_of_iterations=self.config['clustering_niterations'])
        if self.config['reference_pdb']:
            self.calculate_rmsd(reference_pdb=self.config['reference_pdb'])
        self.save_config()
        self.draw_plots()
        self.save_models(replicas=self.config['save_replicas'], topn=self.config['save_topn'],
                         clusters=self.config['save_clusters'], medoids=self.config['save_medoids'])

    def save_config(self):
        if self.config['save_config_file']:
            with open(self.config['work_dir']+'/config.ini', 'w') as configfile:
                for k in self.config:
                    peptide_counter = 0
                    value = self.config[k]
                    line = str(k)+': '
                    if k == 'ligand':
                        for lgnd in value:
                            if peptide_counter==0:
                                line = 'peptide: '
                            else:
                                line += '\nadd-peptide: '
                            for element in lgnd:
                                line += str(element)+' '
                            peptide_counter += 1
                    elif isinstance(value, tuple):
                        line = str(k) + ': '
                        for item in value:
                            line += str(item)+' '
                    else:
                        line = str(k) + ': '+str(value)
                    configfile.write('\n'+line)

    def setup_job(self):
        print('CABS-docking job {0}'.format(self.config['receptor']))
        # Preparing the initial complex
        print(' Building complex...')
        self.initial_complex = ProteinComplex(self.config)
        print(' ... done.')

    def setup_cabs_run(self):
        # Initializing CabsRun instance
        self.cabsrun = CabsRun(self.initial_complex, self.prepare_restraints(), self.config)
        return self.cabsrun

    def execute_cabs_run(self):
        print('CABS simulation starts.')
        self.cabsrun.run()
        print('CABS simuation is DONE.')

    def load_output(self, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param ftraf: path to TRAF file
        :param fseq: path to SEQ file
        :return: returns trajectory.Trajectory instance
        """
        print("load_output")
        if ftraf is not None and fseq is not None:
            self.trajectory = Trajectory.read_trajectory(ftraf, fseq)
        else:
            self.trajectory = self.cabsrun.get_trajectory()
        self.trajectory.number_of_peptides = len(self.config['ligand'])
        self.trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
        self.trajectory.align_to(self.initial_complex.receptor)
        return self.trajectory

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        print("score_results")
        # Filtering the trajectory
        self.filtered_trajectory, self.filtered_ndx = Filter(self.trajectory, n_filtered).cabs_filter()
        # Clustering the trajectory
        self.medoids, self.clusters_dict, self.clusters = Clustering(
            self.filtered_trajectory,
            'chain ' + ','.join(
                self.initial_complex.ligand_chains,
            )
        ).cabs_clustering(number_of_medoids=number_of_medoids, number_of_iterations=number_of_iterations)

    def calculate_rmsd(self, reference_pdb=None, save=True):
        print('calculate_rmsd')
        if save:
            odir = self.config['work_dir'] + '/output_data'
            try:
                mkdir(odir)
            except OSError:
                pass
        all_results = {}
        for pept_chain in self.initial_complex.ligand_chains:
            aln_path = None if not save else self.config[
                                                 'work_dir'] + '/output_data/target_alignment_%s.csv' % pept_chain
            self.rmslst[pept_chain] = self.trajectory.rmsd_to_reference(
                self.initial_complex.receptor_chains,
                ref_pdb=reference_pdb,
                pept_chain=pept_chain,
                align_mth=self.config['align'],
                alignment=self.config['reference_alignment'],
                path=aln_path,
                pept_align_kwargs={'fname': self.config['reference_alignment']},
                target_align_kwargs={'fname': self.config['reference_alignment']}
            )
            rmsds = [header.rmsd for header in self.medoids.headers]
            results = {}
            results['rmsds_all'] = [header.rmsd for header in self.trajectory.headers]
            results['rmsds_filtered'] = [header.rmsd for header in self.filtered_trajectory.headers]
            results['rmsds_medoids'] = rmsds
            results['lowest_all'] = sorted(results['rmsds_all'])[0]
            results['lowest_filtered'] = sorted(results['rmsds_filtered'])[0]
            results['lowest_medoids'] = sorted(results['rmsds_medoids'])[0]
            # Saving rmsd results
            if save:
                with open(odir + '/lowest_rmsds_%s.txt' % pept_chain, 'w') as outfile:
                    outfile.write(
                        'lowest_all; lowest_filtered; lowest_medoids\n {0};{1};{2}'.format(results['lowest_all'],
                                                                                           results['lowest_filtered'],
                                                                                           results['lowest_medoids'], )
                    )
                for type in ['all', 'filtered', 'medoids']:
                    with open(odir + '/{0}_rmsds_{1}.txt'.format(type, pept_chain), 'w') as outfile:
                        for rmsd in results['rmsds_' + type]:
                            outfile.write(str(rmsd) + ';\n')
            all_results[pept_chain] = results
        return all_results

    def draw_plots(self, plots_dir=None):
        print('draw_plots')
        # set the plots dir
        if plots_dir is None:
            pltdir = self.config['work_dir'] + '/plots'
            try:
                mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir

        graph_RMSF(self.trajectory, self.initial_complex.receptor_chains, pltdir + '/RMSF')

        # RMSD-based graphs
        if self.config['reference_pdb']:
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD([self.trajectory, self.filtered_trajectory],
                            [rmslst, rmslst[self.filtered_ndx,]],
                            ['all models', 'top 1000 models'],
                            pltdir + '/E_RMSD_%s' % k)
                plot_RMSD_N(rmslst.reshape(self.config['replicas'], -1),
                            pltdir + '/RMSD_frame_%s' % k)

        # Contact maps
        if self.config['contact_maps']:
            self.mk_cmaps(self.trajectory, self.medoids, self.clusters_dict, self.filtered_ndx, 4.5, pltdir)

    def save_models(self, replicas=True, topn=True, clusters=True, medoids='AA'):
        # output folder
        output_folder = self.config['work_dir'] + '/output_pdbs'
        try:
            mkdir(output_folder)
        except OSError:
            pass
        print('save_models')
        # Saving the trajectory to PDBs:
        if replicas:
            self.trajectory.to_pdb(mode='replicas', to_dir=output_folder)
        # Saving top1000 models to PDB:
        if topn:
            self.filtered_trajectory.to_pdb(mode='replicas', to_dir=output_folder, name='top1000')
        # Saving clusters in CA representation
        if clusters:
            for i, cluster in enumerate(self.clusters):
                cluster.to_pdb(mode='replicas', to_dir=output_folder, name='cluster_{0}'.format(i))
        # Saving top10 models:
        if medoids == 'CA':
            # Saving top 10 models in CA representation:
            self.medoids.to_pdb(mode='models', to_dir=output_folder, name='model')
        elif medoids == 'AA':
            # Saving top 10 models in AA representation:
            pdb_medoids = self.medoids.to_pdb()
            if self.config['AA_rebuild']:
                from cabsDock.ca2all import ca2all
                for i, fname in enumerate(pdb_medoids):
                    ca2all(fname, output=output_folder + '/' + 'model_{0}.pdb'.format(i), iterations=1,
                           verbose=False)

    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, plots_dir):
        scmodeler = SCModeler(ca_traj.template)
        sc_traj_full = scmodeler.calculate_sc_traj(ca_traj.coordinates)
        sc_traj_1k = sc_traj_full.reshape(1, -1, len(ca_traj.template), 3)[:, top1k_inds, :, :]
        sc_med = scmodeler.calculate_sc_traj(meds.coordinates)
        shp = sc_med.shape
        sc_med = sc_med.reshape((shp[1], shp[0]) + shp[2:])

        cmapdir = self.config['work_dir'] + '/contact_maps'
        try:
            mkdir(cmapdir)
        except OSError:
            pass
        rchs = self.initial_complex.receptor_chains
        lchs = self.initial_complex.ligand_chains

        targ_cmf = ContactMapFactory(rchs, rchs, ca_traj.template)
        cmfs = {lig: ContactMapFactory(rchs, lig, ca_traj.template) for lig in lchs}
        cmap10ktarg = reduce(operator.add, targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg.zero_diagonal()
        cmap10ktarg.save_all(cmapdir + '/target_all', break_long_x=0, norm_n=True)

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

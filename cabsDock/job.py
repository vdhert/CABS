"""
Module for running cabsDock jobs.
"""

import re
import operator
from os import getcwd, mkdir
from os.path import exists, isdir, join, abspath

from cabsDock.cluster import Clustering
from protein import ProteinComplex
from restraints import Restraints
from cabs import CabsRun
from utils import ProgressBar
from cabsDock.utils import SCModeler
from cabsDock.utils import _chunk_lst
from cabsDock.plots import plot_E_RMSD
from cabsDock.plots import plot_RMSD_N
from cabsDock.plots import graph_RMSF
from utils import check_peptide_sequence
from trajectory import Trajectory
from cabsDock.cmap import ContactMapFactory
from filter import Filter

__all__ = ['Job']


class Config(dict):
    """
    Smart dictionary that can append items with 'ligand' key instead of overwriting them.
    TODO: universal list of associative options assigned to specific keywords: ligand, restraints etc.
    """
    def __init__(self, config):
        dict.__init__(self, config)

    def read_file(self, filename):
        """Config file parser"""
        patt = re.compile('(\w+)\s+=\s+([\w:.-]+)')
        patt_start = re.compile('start_(\w+)', re.IGNORECASE)
        patt_stop = re.compile('stop', re.IGNORECASE)
        patt_lig = re.compile('ligand\s+=\s+(.*)', re.IGNORECASE)
        with open(filename) as f:
            lines = iter([re.sub('\s+', ' ', line.partition('#')[0]) for line in f.readlines()])

        for line in lines:
            match = re.match(patt_start, line)
            if match:
                key = match.group(1).lower()
                self[key] = []
                try:
                    line = next(lines)
                    while not re.match(patt_stop, line):
                        self[key].append(line.strip())
                        line = next(lines)
                except StopIteration:
                    pass
            else:
                match = re.match(patt_lig, line)
                if match:
                    lig = tuple(w.strip() for w in match.group(1).split(','))
                    if 'ligand' not in self:
                        self['ligand'] = []
                    self['ligand'].append(lig)
                else:
                    match = re.match(patt, line)
                    if match:
                        key = match.group(1).lower()
                        if '.' in match.group(2):
                            try:
                                val = float(match.group(2))
                            except ValueError:
                                val = match.group(2)
                        else:
                            try:
                                val = int(match.group(2))
                            except ValueError:
                                if re.match('\A(yes|y|true|t)\Z', match.group(2), re.IGNORECASE):
                                    val = True
                                elif re.match('\A(no|n|false|f)\Z', match.group(2), re.IGNORECASE):
                                    val = False
                                else:
                                    val = match.group(2)
                        self[key] = val

    def add_config(self, config):
        """Parser for job parameters"""
        if 'ligand' in config:
            if 'ligand' not in self:
                self['ligand'] = []
            val = config['ligand']
            if type(val) in (list, tuple):
                self['ligand'].extend(val)
            elif type(val) is str:
                self['ligand'].append((val,))
            del config['ligand']
        self.update(config)
        return self

    def fix_ligands(self):
        """Function that makes sure each ligand entry in config['ligand'] list is a 3-element tuple"""
        if 'ligand' in self:
            for i, ligand in enumerate(self['ligand']):
                for j in range(len(ligand), 3):
                    self['ligand'][i] += ('random',)
        return self


class Job:
    """
    Class representing single cabsDock job.
    """

    def __repr__(self):
        return '\n'.join([k + ' : ' + str(v) for k, v in sorted(self.config.items())])

    def __init__(self, **kwargs):
        """
        Job can be initialized by:
        1. receptor=[receptor_input] The only required parameter, other taken from 'defaults'.
        2. By providing location of the config file as in 'config=[path_to_file]'.
        3. All parameters can be overwritten by specifying parameter=[value].
        """

        defaults = {
            'work_dir': getcwd(),
            'replicas': 10,
            'mc_cycles': 50,
            'mc_steps': 50,
            't_init': 2.0,
            't_final': 1.0,
            'replicas_dtemp': 0.5,
            'initial_separation': 20.0,
            'ligand_insertion_clash': 0.5,
            'ligand_insertion_attempts': 1000,
            'ca_restraints_strength': 1.0,
            'sg_restraints_strength': 1.0,
            'receptor_restraints': (4, 5.0, 15.0),  # sequence gap, min length, max length
            'dssp_command': 'mkdssp',
            'fortran_compiler': ('gfortran', '-O2'),  # build (command, flags)
            'filtering': 1000,  # number of models to filter
            'clustering_nmedoids': 10,
            'clustering_niterations': 100,  # number of clusters, iterations
            'benchmark': False,
            'AA_rebuild': True,
            'contact_maps': True,
            'reference_pdb': None,
            'align': 'SW',
            'reference_alignment': None,
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

        # Config processing:
        self.config = Config(defaults)

        # check if config should be read from a file
        if 'config' in kwargs:
            if exists(kwargs['config']):
                self.config.read_file(kwargs['config'])
            else:
                raise IOError('Config file: ' + kwargs['config'] + ' does not exist!')

        # update config with kwargs
        self.config.add_config(kwargs).fix_ligands()

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

    def cabsdock(self, withcabs=True, ext_old_ids=None, ext_initial_complex=None, ftraf=None, fseq=None):
        if withcabs:
            self.setup()
            self.execute()
            initial_complex = self.initial_complex
            old_ids = initial_complex.receptor.old_ids
        else:
            initial_complex = ext_initial_complex
            old_ids = ext_old_ids
        self.load_output(old_ids, initial_complex, ftraf, fseq)
        self.score_results(n_filtered=self.config['filtering'], number_of_medoids=self.config['clustering_nmedoids'], number_of_iterations=self.config['clustering_niterations'])
        if self.config['reference_pdb']:
            self.calculate_rmsd(reference_pdb=self.config['reference_pdb'])
        self.draw_plots()
        self.save_models()

    def setup(self):
        print('CABS-docking job {0}'.format(self.config['receptor']))
        # Preparing the initial complex
        print(' Building complex...')
        self.initial_complex = ProteinComplex(self.config)
        print(' ... done.')

        # Generating the restraints
        restraints = \
            Restraints(self.initial_complex.receptor.generate_restraints(*self.config['receptor_restraints']))
        add_restraints = Restraints(self.config.get('ca_restraints'))
        add_restraints += Restraints(self.config.get('sg_restraints'), sg=True)
        add_restraints += Restraints(self.config.get('ca_restraints_file'))
        add_restraints += Restraints(self.config.get('ca_restraints_file'), sg=True)
        restraints += add_restraints.update_id(self.initial_complex.new_ids)

        # Initializing CabsRun instance
        self.cabsrun = CabsRun(self.initial_complex, restraints, self.config)
        return self.cabsrun

    def execute(self):
        print('CABS simulation starts.')
        self.cabsrun.run()
        print('CABS simuation is DONE.')

    def load_output(self, old_ids, initial_complex, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param number_of_peptides:
        :param old_ids:
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
        self.trajectory.template.update_ids(old_ids, pedantic=False)
        self.trajectory.align_to(initial_complex.receptor)
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
            aln_path = None if not save else self.config['work_dir'] + '/output_data/target_alignment_%s.csv' % pept_chain
            #~ self.trajectory.rmsd_to_native_test(
                #~ reference_pdb,
                #~ self.initial_complex.receptor_chains,
                #~ pept_chain,
                #~ pept_chain,
                #~ )
            self.rmslst[pept_chain] = self.trajectory.rmsd_to_reference(
                self.initial_complex.receptor_chains,
                ref_pdb=reference_pdb,
                pept_chain=pept_chain,
                align_mth=self.config['align'],
                alignment=self.config['reference_alignment'],
                path=aln_path
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
                with open(odir+'/rmsds_%s.txt' % pept_chain, 'w') as outfile:
                    outfile.write(
                        'lowest_all; lowest_filtered; lowest_medoids\n {0};{1};{2}'.format(results['lowest_all'],
                                                                                           results['lowest_filtered'],
                                                                                           results['lowest_medoids'], )
                    )
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
                             pltdir + '/E_RMSD_%s' % k)
                plot_RMSD_N(rmslst.reshape(self.config['replicas'], -1),
                            pltdir + '/RMSD_frame_%s' % k)

        # Contact maps
        if self.config['contact_maps']:
            self.mk_cmaps(self.trajectory, self.medoids, self.clusters_dict, self.filtered_ndx, 4.5, pltdir)

    def save_models(self, replicas=True, topn=True, clusters=True, medoids='AA'):
        #output folder
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
        cmap10ktarg.save_all(cmapdir + '/target_all', break_long_x=0)

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

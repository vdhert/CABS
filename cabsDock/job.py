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
from utils import kmedoids
from utils import SCModeler
from utils import check_peptide_sequence
from trajectory import Trajectory
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap
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
            if type(val) is list:
                self['ligand'].extend(val)
            elif type(val) is tuple:
                self['ligand'].append(val)
            elif type(val) is str:
                self['ligand'].append((val,))
            del config['ligand']
            # checks the input peptide sequence for non-standard amino acids.
            [check_peptide_sequence(peptide[0]) for peptide in self['ligand']]
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
            'clustering': (10, 100),  # number of clusters, iterations
            'native_pdb': None
        }

        self.config = Config(defaults)

        # check if config should be read from a file
        if 'config' in kwargs:
            if exists(kwargs['config']):
                self.config.read_file(kwargs['config'])
            else:
                raise IOError('Config file: ' + kwargs['config'] + ' does not exist!')

        # update config with kwargs
        self.config.add_config(kwargs).fix_ligands()

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

    def run_job(self):
        work_dir = self.config['work_dir']
        print('CABS-docking job {0}'.format(self.config['receptor']))
        # prepare initial complex
        # noinspection PyAttributeOutsideInit
        print(' Building complex...')
        self.initial_complex = ProteinComplex(self.config)

        print(' ... done.')
        # generate restraints
        # noinspection PyAttributeOutsideInit
        self.restraints = \
            Restraints(self.initial_complex.receptor.generate_restraints(*self.config['receptor_restraints']))
        add_restraints = Restraints(self.config.get('ca_restraints'))
        add_restraints += Restraints(self.config.get('sg_restraints'), sg=True)
        add_restraints += Restraints(self.config.get('ca_restraints_file'))
        add_restraints += Restraints(self.config.get('ca_restraints_file'), sg=True)
        self.restraints += add_restraints.update_id(self.initial_complex.new_ids)

        # run cabs
        print('CABS simulation starts.')
        cabs_run = CabsRun(self.initial_complex, self.restraints, self.config)
        cabs_run.run()
        # bar = ProgressBar(100, msg='CABS is running:')
        # while cabs_run.is_alive():
        #     print('isAlive')
        #     bar.update(cabs_run.status())
        #     sleep(5)
        # bar.done()
        print('CABS simuation is DONE.')
        if self.config['native_pdb']:
            print('Calculating RMSD to the native structure...')
            trajectory = cabs_run.get_trajectory()
            trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
            print(
                'The native complex loaded from {0} consists of receptor (chain(s) {1}) and peptide(s) (chains(s) {2}).'
                .format(
                    self.config['native_pdb'],
                    self.config['native_receptor_chain'],
                    self.config['native_peptide_chain']
                )
            )
            trajectory.rmsd_to_native(native_pdb=self.config['native_pdb'],
                                      native_receptor_chain=self.config['native_receptor_chain'],
                                      native_peptide_chain=self.config['native_peptide_chain'],
                                      model_peptide_chain=self.initial_complex.ligand_chains[0])
            trajectory.align_to(self.initial_complex.receptor)
        else:
            trajectory = cabs_run.get_trajectory()
            trajectory.align_to(self.initial_complex.receptor)
            trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
        tra, flt_inds = Filter(trajectory).filter()
        # MC: Functionality moved to a separate class cabsDock.clustering.Clustering (IN PROGRESS)
        medoids, clusters = Clustering(tra, 'chain ' + ','.join(self.initial_complex.ligand_chains)).cabs_clustering()

        #~ import pickle
        #~ with open("traj.pck", "w") as f:
            #~ pickle.dump(trajectory, f)
        #~ with open("clst.pck", "w") as f:
            #~ pickle.dump(clusters, f)
        #~ with open("flti.pck", "w") as f:
            #~ pickle.dump(flt_inds, f)

        self.mk_cmaps(trajectory, clusters, flt_inds, 4.5)

        #Saving the models to PDB
        # for i, medoid in enumerate(medoids.coordinates[0]):
        #     filename = join(work_dir, 'model_%d.pdb' % i)
        #
        #
        #
        #
        #     m.save_to_pdb(filename, bar_msg='Saving %s' % filename)
        #
        # for i, m in enumerate(trajectory.coordinates, 1):
        #     filename = join(work_dir, 'replica_%d.pdb' % i)
        #     replica = Trajectory(trajectory.template, m, None).to_atoms()
        #     replica.save_to_pdb(filename, bar_msg='Saving %s' % filename)

        # dictionary holding results to be returned for use in the Benchmark class
        rmsds = [header.rmsd for header in medoids.headers ]
        results = {}
        results['rmsds_10k'] = [header.rmsd for header in trajectory.headers]
        results['rmsds_1k'] = [header.rmsd for header in tra.headers]
        results['rmsds_10'] = rmsds
        results['lowest_10k'] = sorted(results['rmsds_10k'])[0]
        results['lowest_1k'] = sorted(results['rmsds_1k'])[0]
        results['lowest_10'] = sorted(results['rmsds_10'])[0]
        print('... done.')
        return results

    def mk_cmaps(self, ca_traj, clusts, top1k_inds, thr):
        scmodeler = SCModeler(self.initial_complex)
        sc_traj_full = scmodeler.calculate_sc_traj(ca_traj.coordinates)

        #~ import imp
        #~ pdbx = imp.load_source('test', '/usr/lib/python2.7/pdb.py')
        #~ pdbx.set_trace()

        cmapdir = self.config['work_dir'] + '/contact_maps'
        try: mkdir(cmapdir)
        except OSError: pass
        rchs = self.initial_complex.receptor_chains
        lchs = self.initial_complex.ligand_chains
        cmfs = {lig: ContactMapFactory(rchs, lig, ca_traj.template) for lig in lchs}
        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(sc_traj_full, thr)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(cmapdir + '/replica_%i_ch_%s' % (n + 1, lig))
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + '/all_ch_%s' % lig)
            sc_traj_1k = sc_traj_full.reshape(1, -1, len(ca_traj.template), 3)[:,top1k_inds,:,:]
            cmap1k = cmf.mk_cmap(sc_traj_1k, thr)[0]
            cmap1k.save_all(cmapdir + '/top1000_ch_%s' % lig)
            for cn, clust in clusts.items():
                ccmap = cmf.mk_cmap(sc_traj_1k, thr, frames=clust)[0]
                ccmap.save_all(cmapdir + '/cluster_%i_ch_%s' % (cn, lig))



if __name__ == '__main__':
    j = Job(receptor='1jbu:H', ligand = [['EEWEVLCWTWETCER']], mc_cycles=10, mc_steps=1, replicas=2, native_pdb='1jbu',
                               native_receptor_chain='H',
                               native_peptide_chain='X')
    print j.run_job()
    # j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=2, mc_steps=2, replicas=2, )
    # j.run_job()

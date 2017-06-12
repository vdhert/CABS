"""
Module for running cabsDock jobs.
"""

import re
import operator
from os import getcwd, mkdir
from os.path import exists, isdir, join, abspath
from time import sleep

import pickle

from protein import ProteinComplex
from restraints import Restraints
from cabs import CabsRun
from utils import ProgressBar, kmedoids
from trajectory import Trajectory
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap

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
            'dssp_command': 'dssp',
            'fortran_compiler': ('gfortran', '-O2'),    # build (command, flags)
            'filtering': 1000,  # number of models to filter
            'clustering': (10, 100)  # number of clusters, iterations
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

        # prepare initial complex
        self.initial_complex = ProteinComplex(self.config)

        if 'dbg' in kwargs: #ROR
            with open("test_complex.pck", "w") as f:
                pickle.dump(self.initial_complex, f)

        # generate restraints
        self.restraints = \
            Restraints(self.initial_complex.receptor.generate_restraints(*self.config['receptor_restraints']))
        add_restraints = Restraints(self.config.get('ca_restraints'))
        add_restraints += Restraints(self.config.get('sg_restraints'), sg=True)
        add_restraints += Restraints(self.config.get('ca_restraints_file'))
        add_restraints += Restraints(self.config.get('ca_restraints_file'), sg=True)
        self.restraints += add_restraints.update_id(self.initial_complex.new_ids)

        # run cabs
        cabs_run = CabsRun(self.initial_complex, self.restraints, self.config)
        cabs_run.start()
        bar = ProgressBar(100, msg='CABS is running:')
        while cabs_run.is_alive():
            bar.update(cabs_run.status())
            sleep(0.1)
        bar.done()
        trajectory = cabs_run.get_trajectory()
        if 'dbg' in kwargs:     #ROR
            with open("test_traj.pck", "w") as f:
                pickle.dump(trajectory, f)
        trajectory.align_to(self.initial_complex.receptor)
        trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
        tra = trajectory.filter(self.config['filtering'])

        #  od tego miejsca poprawic
        lig_chains = ','.join(self.initial_complex.ligand_chains)
        if lig_chains:
            ligs = tra.select('chain %s' % lig_chains)
        else:
            ligs = tra
        D = ligs.rmsd_matrix(msg='Calculating rmsd matrix')
        M, C = kmedoids(D, *self.config['clustering'])

        #TO-start: cmap factory init; cmaps for replicas
        cmapdir = self.config['work_dir'] + '/contact_maps'
        try:
            mkdir(cmapdir)
        except OSError:
            pass
        cmfs = {lig: ContactMapFactory(self.initial_complex.receptor_chains, lig, trajectory.template) for lig in self.initial_complex.ligand_chains}
        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(trajectory.coordinates, 6.5)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(cmapdir + '/replica_%i_ch_%s' % (n + 1, lig))
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + '/all_ch_%s' % lig)
            cmap1k = cmf.mk_cmap(tra.coordinates, 6.5)[0]
            cmap1k.save_all(cmapdir + '/top1000_ch_%s' % lig)
            for cn, clust in C.items():
                ccmap = cmf.mk_cmap(tra.coordinates, 6.5, frames=clust)[0]
                ccmap.save_all(cmapdir + '/cluster_%i_ch_%s' % (cn, lig))
        #TO-end

        if 'dbg' in kwargs:     #ROR
            with open("test_clusters.pck", "w") as f:
                pickle.dump(C, f)
        medoids = [tra.get_model(m) for m in M]
        for i, m in enumerate(medoids, 1):
            filename = join(work_dir, 'model_%d.pdb' % i)
            m.save_to_pdb(filename, bar_msg='Saving %s' % filename)

        for i, m in enumerate(trajectory.coordinates, 1):
            filename = join(work_dir, 'replica_%d.pdb' % i)
            replica = Trajectory(trajectory.template, m, None).to_atoms()
            replica.save_to_pdb(filename, bar_msg='Saving %s' % filename)


if __name__ == '__main__':
    j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=50,  mc_steps=1, replicas=10)

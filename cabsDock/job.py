"""
Module for running cabsDock jobs.
"""

import re
from os import getcwd, mkdir
from os.path import exists, isdir
from time import sleep

from protein import ProteinComplex
from restraints import Restraints
from cabs import CabsRun
from utils import ProgressBar

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
        trajectory.align_to(self.initial_complex.receptor)
        trajectory.template.update_ids(self.initial_complex.receptor.old_ids, pedantic=False)
        tra = trajectory.filter(27)
        tra.to_atoms().save_to_pdb('dupa.pdb', bar_msg='Saving ...')


if __name__ == '__main__':
    j = Job(receptor='1rjk:A', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=3,  mc_steps=3, replicas=1)


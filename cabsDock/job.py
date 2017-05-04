"""
Module for running cabsDock jobs.
"""

from os import getcwd, mkdir
from os.path import exists, isdir
import re
from protein import *


class Job:
    """
    Class representing single cabsDock job.
    """

    class Config(dict):
        """Class handling job configuration"""
        def __init__(self):
            dict.__init__(self)

        def update(self, other):
            dict.update(self, other)

    @staticmethod
    def read_config(config_file):
        """Reads and parses config file. Returns a dictionary with options."""
        config = {}
        patt = re.compile('(\w+)\s+=\s+([\w:.-]+)')
        patt_start = re.compile('start_(\w+)', re.IGNORECASE)
        patt_stop = re.compile('stop', re.IGNORECASE)
        patt_lig = re.compile('ligand\s+=\s+(.*)', re.IGNORECASE)

        with open(config_file) as f:
            lines = iter([re.sub('\s+', ' ', line.partition('#')[0]) for line in f.readlines()])

        for line in lines:
            match = re.match(patt_start, line)
            if match:
                key = match.group(1).lower()
                config[key] = []
                try:
                    line = next(lines)
                    while not re.match(patt_stop, line):
                        config[key].append(line.strip())
                        line = next(lines)
                except StopIteration:
                    pass
            else:
                match = re.match(patt_lig, line)
                if match:
                    if 'ligand' not in config:
                        config['ligand'] = []
                    config['ligand'].append([w.strip() for w in match.group(1).split(',')])
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
                        config[key] = val
        return config

    def __repr__(self):
        return '\n'.join([k + ' : ' + str(v) for k, v in sorted(self.config.items())])

    def __init__(self, **kwargs):
        """
        Job can be initialized by:
        1. None - therefore taking default parameters stored in self.config,
           or if config.txt exist in current working directory it overwrites default values.
        2. By providing location of the config file as config=[filepath], or providing path
           to the directory containing config.txt file as work_dir=[dirpath].
        3. All parameters can be overwritten by specifying parameter=[value].
        """

        self.config = Config()

        default_config = {
            'work_dir': getcwd(),
            'replicas': 10,
            'mc_cycles': 50,
            'mc_steps': 50,
            't_init': 2.0,
            't_final': 1.0,
            'replicas_dtemp': 0.5,
            'initial_separation': 20.0,
            'ligand_insertion_clash': 2.5,
            'ligand_insertion_attempts': 1000,
            'ca_restraints_strength': 1.0,
            'sg_restraints_strength': 1.0
        }

        # reading config file
        if 'config' in kwargs:
            if exists(kwargs['config']):
                config_from_file = self.read_config(kwargs['config'])
            else:
                raise IOError('Config file: ' + kwargs['config'] + ' does not exist!')
        self.config.update(kwargs)

        # checks if work_dir exists, creates it otherwise
        work_dir = self.config['work_dir']
        if exists(work_dir):
            if not isdir(work_dir):
                raise Exception('File %s already exists and is not a directory' % work_dir)
        else:
            mkdir(work_dir)

        self.initial_complex = ProteinComplex(self.config)


if __name__ == '__main__':
    j = Job(receptor='1rjk:A', config='../test/config.txt')
    print j
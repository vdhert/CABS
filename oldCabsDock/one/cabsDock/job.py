"""
Module for running cabsDock jobs.
"""

from os import getcwd, mkdir
from os.path import exists, isdir, join
from urllib2 import HTTPError
from pdb import *
from atom import *
from protein import *
from cabs import *
import re


class Job:
    """
    Class representing single cabsDock job.
    """

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

    LATTICE = CabsLattice()  # define a lattice object for casting

    def __init__(self, **kwargs):
        """
        Job can be initialized by:
        1. None - therefore taking default parameters stored in self.config,
           or if config.txt exist in current working directory it overwrites default values.
        2. By providing location of the config file as config=[filepath], or providing path
           to the directory containing config.txt file as work_dir=[dirpath].
        3. All parameters can be overwritten by specifying parameter=[value].
        """

        self.config = {
            'work_dir': getcwd(),
            'receptor_input': 'input.pdb',
            'dssp_command': 'dssp',
            'replicas': 10,
            'mc_cycles': 50,
            'mc_steps': 50,
            't_init': 2.0,
            't_final': 1.0,
            'replicas_dtemp': 0.5,
            'brake_chains_at': 4.5
        }

        # reading config file
        if 'config' in kwargs:
            if exists(kwargs['config']):
                self.config.update(self.read_config(kwargs['config']))
            else:
                raise IOError('Config file: ' + kwargs['config'] + ' does not exist!')
        self.config.update(kwargs)

        # checks if work_dir exists, creates it otherwise
        self.work_dir = self.config['work_dir']
        if exists(self.work_dir):
            if not isdir(self.work_dir):
                raise IOError('File ' + self.work_dir + ' already exists and is not a directory')
        else:
            mkdir(self.work_dir)

        receptor = self.prepare_receptor(self.config['receptor'])
        self.original_ids = receptor.fix_broken_chains(cut_off=self.config['brake_chains_at'])
        self.new_ids = {v: k for k, v in self.original_ids.iteritems()}
        self.protein_complex = Complex(receptor, self.config['replicas'])

        if 'ligand' in self.config:
            for lig in self.config['ligand']:
                for i in range(len(lig), 3):
                    lig.append('random')
                ligand, rnd = self.prepare_ligand(lig[0])
                if lig[1] == 'random':
                    rnd = True
                loc = self.fix_initial_location(lig[2])
                self.protein_complex.add_ligand(ligand, randomize=rnd, location=loc)

        # ca_restraints = self.prepare_ca_restraints()
        # sg_restraints = self.prepare_sg_restraints()

    def prepare_receptor(self, receptor_input):
        """
        Returns a Protein object with single copy of the receptor in CA-only representation
        depending on the content of the [receptor_input] parameter.
        """
        # checks if receptor_input is a file and if it exists, loads it into a Pdb object
        if exists(receptor_input):
            receptor_pdb = Pdb(pdb_file=receptor_input)
        elif exists(join(self.work_dir, receptor_input)):
            receptor_pdb = Pdb(pdb_file=join(self.work_dir, receptor_input))

        # or if receptor_input is a valid PDB code and downloads it
        else:
            try:
                receptor_pdb = Pdb(pdb_code=receptor_input[:4])
            except HTTPError:
                raise IOError(receptor_input + ' is not a valid PDB Code')
            m = re.match('.{4}:([A-Z]*)', receptor_input)
            if m:
                select_chains = m.group(1)

        selection = 'name CA and not HETERO'
        if 'select_chains' in locals():
            selection += ' and chain ' + ','.join(select_chains)
        rec_ca = Atoms(receptor_pdb).remove_alternative_locations().select(selection).models()[0]
        rec_ca.update_sec(receptor_pdb.dssp(self.config['dssp_command']))
        return Protein(rec_ca)

    def prepare_ligand(self, ligand_input):
        """
        Returns a tuple(Protein, bool).
        Protein contains ligand coordinates in CA-only representation read from file,
        downloaded from the PDB, or randomly generated, depending of the content of the [ligand_input].
        bool is True if conformation should be randomized, False otherwise.
        """
        # checks if ligand_input is a file and if it exists, loads it into a Pdb object
        if exists(ligand_input):
            ligand_pdb = Pdb(pdb_file=ligand_input)
        elif exists(join(self.work_dir, ligand_input)):
            ligand_pdb = Pdb(pdb_file=join(self.work_dir, ligand_input))

        # or if ligand_pdb is a valid PDB code and downloads it
        else:
            try:
                ligand_pdb = Pdb(pdb_code=ligand_input[:4])
                if ':' in ligand_input:
                    select_chain = ligand_input.split(':')[1]

            # if that fails ligand_input is considered a SEQUENCE[:SECONDARY]
            except HTTPError:
                if ':' in ligand_input:
                    seq, sec = ligand_input.upper().split(':')
                else:
                    seq = ligand_input.upper()

        if 'ligand_pdb' in locals():
            selection = 'name CA and not HETERO'
            if 'select_chain' in locals():
                selection += ' and chain ' + select_chain
            return Protein(Atoms(ligand_pdb).remove_alternative_locations().select(selection)), False

        else:
            if not re.search('[^A-Z]', seq):
                if 'sec' in locals():
                    return Protein(seq + ':' + sec).random_conformation(), True
                else:
                    return Protein(seq).random_conformation(), True
            else:
                raise ValueError('Invalid ligand input: ' + ligand_input)

    def fix_initial_location(self, loc):
        if loc == 'random' or loc == 'keep':
            return loc
        else:
            return '+'.join([self.new_ids[l] for l in loc.split('+')])

if __name__ == '__main__':
    j = Job(config='../test/config.txt', work_dir='../test', temp=3)
    for k, v in sorted(j.config.iteritems()):
        print k, ":", v
    with open('dupa.pdb', 'w') as f:
        f.write(j.protein_complex.make_pdb())
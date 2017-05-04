from string import ascii_uppercase
import numpy as np
import re
from os.path import abspath, join


__all__ = [
    'CABS_HOME',
    'CABS_SS',
    'RANDOM_LIGAND_LIBRARY',
    'AA_NAMES',
    'aa_to_long',
    'aa_to_short',
    'next_letter'
]

# Location of the CABSdock dir
CABS_HOME = abspath(join(__file__, '../..'))

# Dictionary for conversion of secondary structure from DSSP to CABS
CABS_SS = {'C': 1, 'H': 2, 'T': 3, 'E': 4, 'c': 1, 'h': 2, 't': 3, 'e': 4}

# Library of 1000 random peptides with up to 50 amino acids each
RANDOM_LIGAND_LIBRARY = np.reshape(np.fromfile(CABS_HOME + '/data/libLig50.dat', sep=' '), (1000, 50, 3))

# Dictionary for amino acid name conversion
AA_NAMES = {
    'A': 'ALA', 'B': 'ASX', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'J': 'XLE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'O': 'HOH',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR',
    'U': 'UNK', 'V': 'VAL', 'W': 'TRP', 'X': 'XAA', 'Y': 'TYR',
    'Z': 'GLX'
}


def aa_to_long(seq):
    """Converts short amino acid name to long."""
    s = seq.upper()
    if s in AA_NAMES:
        return AA_NAMES[s]
    else:
        raise InvalidAAName(seq, 1)


def aa_to_short(seq):
    """Converts long amino acid name to short."""
    s = seq.upper()
    for short, full in AA_NAMES.items():
        if full == s:
            return short
    else:
        raise InvalidAAName(seq, 3)


def next_letter(taken_letters):
    """Returns next available letter for new protein chain."""
    return re.sub('[' + taken_letters + ']', '', ascii_uppercase)[0]


class InvalidAAName(Exception):
    """Exception raised when invalid amino acid name is used"""
    def __init__(self, name, l):
        self.name = (name, l)

    def __str__(self):
        return '%s is not a valid %d-letter amino acid code' % self.name

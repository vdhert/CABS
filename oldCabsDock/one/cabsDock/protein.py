"""
Module for handling proteins and peptides in CA representation.
"""

import re
from atom import *
from vector3d import Vector3d
from string import ascii_uppercase
from copy import deepcopy as cp
from random import randint, random
from pdb import *
from math import sin, cos, pi


__all__ = ['Protein', 'Complex']


class Protein(Atoms):
    """
    Class for proteins. Inherits from Atoms. Has some protein specific methods and attributes.
    """
    AA_NAMES = {
        'A': 'ALA', 'B': 'ASX', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'J': 'XLE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'O': 'HOH',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR',
        'U': 'UNK', 'V': 'VAL', 'W': 'TRP', 'X': 'XAA', 'Y': 'TYR',
        'Z': 'GLX'
    }

    CABS_SS = {'C': 1, 'H': 2, 'T': 3, 'E': 4}

    @classmethod
    def aa_to_long(cls, seq):
        return cls.AA_NAMES[seq]

    @classmethod
    def aa_to_short(cls, seq):
        for short, full in cls.AA_NAMES.items():
            if full == seq:
                return short

    @staticmethod
    def read_peptide_library(filename, pep_length):
        lib = []
        pep = []
        for line in open(filename, 'r').readlines():
            pep.append(Vector3d(line))
            if len(pep) == pep_length:
                lib.append(cp(pep))
                pep = []
        return lib

    MAX_PEPTIDE_LENGTH = 50
    PEPTIDE_LIBRARY = read_peptide_library.__func__('../data/libLig50.dat', MAX_PEPTIDE_LENGTH)

    @classmethod
    def get_random_coordinates(cls, pep_length):
        if pep_length > cls.MAX_PEPTIDE_LENGTH:
            raise ValueError('Max peptide length is %i aa' % cls.MAX_PEPTIDE_LENGTH)
        rand_model = randint(0, len(cls.PEPTIDE_LIBRARY) - 1)
        rand_index = randint(0, cls.MAX_PEPTIDE_LENGTH - pep_length)
        return cp(cls.PEPTIDE_LIBRARY[rand_model][rand_index:rand_index + pep_length])

    def __init__(self, arg=None):
        """
        Initialize with:
        1. SEQUENCE[:SECONDARY]
        2. int for polyalanine
        3. Atoms()
        """
        if type(arg) is str:
            self.atoms = []
            if ':' in arg:
                seq, sec = arg.split(':')
                if len(sec) != len(seq):
                    sec = 'C' * len(seq)
            else:
                seq = arg
                sec = 'C' * len(seq)

            for i, ch in enumerate(seq):
                self.atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i+1,
                        name=' CA ',
                        resname=self.aa_to_long(ch),
                        resnum=i+1,
                        occ=self.CABS_SS[sec[i]]
                    )
                )
        elif type(arg) is int:
            self.atoms = []
            for i in range(arg):
                self.atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i+1,
                        name=' CA ',
                        resname='ALA',
                        resnum=i+1
                    )
                )
        else:
            Atoms.__init__(self, arg)

    def fix_broken_chains(self, cut_off=4.5, used_letters=""):
        """
        Checks for gaps in protein chains (Ca-Ca distance > cut_off). Splits broken chains
        on gaps taking next available letter for the new chain, except for those in used_letters.
        Returns a dictionary with residue ids (new -> old).
        """

        used_letters += ''.join(self.list_chains().keys())

        old_ids = {}
        for chain in self.chains():
            prev = None
            for residue in chain.residues():
                ca = residue.select('name CA')[0]
                res_id = ca.resid_id()
                if prev:
                    chid = prev.chid
                    d = (ca.coord - prev.coord).length()
                    if d > cut_off:
                        chid = sorted(re.sub('[' + used_letters + ']', '', ascii_uppercase))[0]
                        used_letters += chid
                    for a in residue:
                        a.chid = chid
                old_ids[ca.resid_id()] = res_id
                prev = ca
        return old_ids

    def random_conformation(self):
        """
        Generates a random conformation for a peptide taken from a library of random peptides.
        """
        if len(self) <= self.MAX_PEPTIDE_LENGTH:
            coord = self.get_random_coordinates(len(self))
            for i, atom in enumerate(self):
                atom.coord = coord[i]
        return self


def next_letter(taken_letters):
    """
    Returns next available letter for new protein chain.
    """
    return re.sub('[' + taken_letters + ']', '', ascii_uppercase)[0]


class Complex(Atoms):
    """
    Class representing protein complex with multiple models.
    """
    CLASH_CUTOFF = 0.5  # minimum distance between any two atoms when adding molecules to the complex.
    SEPARATION = 20.0   # distance from approximate protein surface at which new molecules are placed
    ATTEMPTS = 100      # maximum number of attempts when adding new ligand in random location

    def __init__(self, receptor, replicas=1):
        """
        Initialized with Atoms object(the receptor). Creates [replicas] copies.
        """
        Atoms.__init__(self)
        self.center_of_mass = receptor.cent_of_mass()
        self.dimension = receptor.max_dimension()
        for i in range(replicas):
            r = cp(receptor)
            self.atoms.extend(r.set_model_number(i + 1))

    def add_ligand(self, ligand, location='random', randomize=False):
        """
        Adds new ligand to all replicas considering other ligands
        to avoid clashes and ambiguity in chain IDs.

        [ligand] is Atoms object with coordinates.
        [location] is a string: either 'random', 'keep' or receptor's assumed binding pocket residues' ids
        joined with '+' i.e. '123:A+127:A+11:B'
        [randomize] is bool telling if coordinates in [ligand] should be randomized
        """
        taken_chains = ''.join(self.list_chains().keys())
        ligand_chid = ligand.atoms[0].chid
        if ligand_chid in taken_chains or ligand_chid not in ascii_uppercase:
            ligand.change_chid(ligand_chid, next_letter(taken_chains))
        loc_vector = None
        lig_com = ligand.cent_of_mass()

        for model in self.models():
            lig = cp(ligand)
            lig.set_model_number(model.atoms[0].model)
            if randomize:
                lig = Protein(lig).random_conformation()
            if location == 'keep':
                lig.move_to(lig_com)
                if model.min_distance(lig) > self.CLASH_CUTOFF:
                    self.atoms.extend(lig)
                else:
                    raise ValueError('Clash!')
            elif location == 'random':
                for i in range(self.ATTEMPTS):
                    theta = random() * pi
                    phi = random() * 2.0 * pi
                    sin_theta = sin(theta)
                    r = self.dimension * 0.5 + self.SEPARATION
                    x = r * sin_theta * cos(phi)
                    y = r * sin_theta * sin(phi)
                    z = r * cos(theta)
                    lig.center_at_origin().move(self.center_of_mass + Vector3d(x, y, z))
                    if model.min_distance(lig) > self.CLASH_CUTOFF:
                        self.atoms.extend(lig)
                        break
                else:
                    raise ValueError('Maximum number of attempts to add new ligand reached.')
            else:
                if not loc_vector:
                    chains = {}
                    for res in location.split('+'):
                        num, chid = res.split(':')
                        if chid in chains:
                            chains[chid].append(num)
                        else:
                            chains[chid] = [num]
                    s = " or ".join(["(chain " + ch + " and resnum " + ",".join(chains[ch]) + ")" for ch in chains])
                    patch = model.select(s)
                    loc_vector = (patch.cent_of_mass() - self.center_of_mass).norm()\
                        * (self.SEPARATION + self.dimension * 0.5)
                lig.center_at_origin().move(self.center_of_mass + loc_vector)
                if model.min_distance(lig) > self.CLASH_CUTOFF:
                    self.atoms.extend(lig)
                else:
                    raise ValueError('Clash!')


if __name__ == '__main__':
    a = Atoms(Pdb(pdb_code='2gb1')).select('name CA')
    c = Complex(a, 10)
    l = Atoms(Pdb(pdb_code='1rjk')).select('chain C and name CA')
    l.change_chid('C', 'A')
    c.add_ligand(l, '16:A+33:A+29:A', True)
    print c.make_pdb()

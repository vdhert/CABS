"""
Module contains two classes: Atom and Atoms. Atom represents single atom/line from the pdb file.
Atoms is a container for Atom objects, without actually specyfing if they are in the same molecule/chain/residue etc.
"""

import re
import numpy as np
from math import sqrt
from copy import deepcopy
from itertools import combinations
from string import ascii_uppercase

from utils import CABS_SS, aa_to_long, smart_flatten, kabsch, ProgressBar
from vector3d import Vector3d


class Atom:
    """
    Class for representation of a single atom.
    """
    def __init__(self, line=None, model=0, **kwargs):
        """
        Constructor. Creates an Atom object from string - ATOM/HETATM line from the pdb file.
        If line is empty creates an empty atom equivalent to:
        Atom('HETATM    0 XXXX XXX X   0       0.000   0.000   0.000  0.00  0.00')
        Passing attribute=value to the constructor overwrites default/read values.
        """
        if line:
            self.model = model
            self.hetatm = (line[:6] == "HETATM")
            self.serial = int(line[6:11])
            self.name = line[11:16].strip()
            self.alt = line[16]
            self.resname = line[17:21].strip()
            self.chid = line[21]
            self.resnum = int(line[22:26])
            self.icode = line[26]
            self.coord = Vector3d(line[27:54])
            self.occ = float(line[54:60])
            self.bfac = float(line[60:66])
            self.tail = line[66:].replace('\n', '')
        else:
            self.model = model
            self.hetatm = True
            self.serial = 0
            self.name = "XXXX"
            self.alt = ""
            self.resname = "XXX"
            self.chid = "X"
            self.resnum = 0
            self.icode = ""
            self.coord = Vector3d()
            self.occ = 0.0
            self.bfac = 0.0
            self.tail = ""

        for arg in kwargs:
            if arg in self.__dict__:
                self.__dict__[arg] = kwargs[arg]

    def __repr__(self):
        line = "ATOM  "
        if self.hetatm:
            line = "HETATM"
        fmt_name = " %-3s" % self.name
        if len(self.name) == 4:
            fmt_name = self.name
        line += "%5d %4s%1s%-4s%1s%4d%1s   %24s%6.2f%6.2f%s" % (
                self.serial,
                fmt_name,
                self.alt,
                self.resname,
                self.chid,
                self.resnum,
                self.icode,
                self.coord,
                self.occ,
                self.bfac,
                self.tail
        )
        return line

    def same_model(self, other):
        """
        Returns True if both atoms belong to the same model. False otherwise.
        """
        return self.model == other.model

    def same_chain(self, other):
        """
        Returns True if both atoms belong to the same chain and model. False otherwise.
        """
        return self.same_model(other) and self.chid == other.chid

    def same_residue(self, other):
        """
        Returns True if both atoms belong to the same residue, chain and model. False otherwise.
        """
        return self.same_chain(other) and self.resnum == other.resnum and self.icode == other.icode

    def is_hydrogen(self):
        """
        Returns true if Atom is hydrogen, false otherwise.
        Determined by the first non-digit character in atom's name. If "H" then hydrogen.
        """
        m = re.search("([A-Z])", self.name)
        return m and m.group(0) == "H"

    def dist2(self, other):
        """
        Returns squared distance between two atoms.
        """
        return (self.coord - other.coord).mod2()

    def distance(self, other):
        """
        Returns distance in Angstroms between two atoms.
        """
        return sqrt(self.dist2(other))

    def min_distance(self, other):
        """
        Returns minimal distance between atom and group of atoms.
        """
        return min(self.distance(atom) for atom in other)

    def match_token(self, token):
        """
        Returns True if Atom matches selection token. False otherwise.
        """
        words = token.split()
        if len(words) == 1:
            # Test here "no argument" selection keywords.
            kword = token.upper()
            if kword == "HETERO":
                return self.hetatm
            else:
                raise Exception('Invalid selection keyword: ' + kword)
        elif len(words) > 1:
            # Here test selection keywords with arguments.
            kword = words[0].upper()
            args = ''.join(words[1:]).split(',')
            if kword == "MODEL":
                return self.model in smart_flatten(args)
            elif kword == "CHAIN":
                return self.chid in args
            elif kword == "RESNUM":
                return self.resnum in smart_flatten(args)
            elif kword == "RESNAME":
                return any([re.match("^%s$" % a, self.resname) for a in args])
            elif kword == "NAME":
                return any([re.match("^%s$" % a, self.name) for a in args])
            else:
                raise Exception('Invalid selection keyword: ' + kword)
        else:
            raise Exception('Invalid selection syntax: ' + token)

    def match(self, sele):
        """
         Returns True if Atom matches selection pattern. False otherwise.
        """
        pattern = deepcopy(sele.tokens)
        for i, t in enumerate(pattern):
            if t.upper() not in Selection.JOINTS:
                pattern[i] = str(self.match_token(t))
        return eval(" ".join(pattern))

    def resid_id(self):
        """Returns a string with residue identification i.e. 123:A"""
        return (str(self.resnum) + self.icode).strip() + ":" + self.chid


class Atoms:
    """
    Container for atoms. Has most methods of a list. Also has methods common
    for all multi-atom objects: move, rotate etc.
    """

    def __init__(self, arg=None):
        """
        Constructor. arg should be either:
            - None
            - list[Atom]
            - any object that has attribute named atoms, which is a list of objects named 'Atom'
            - string 'SEQUENCE' or 'SEQUENCE:SECONDARY' [SECONDARY is H/E/C/T for helix/sheet/turn/coil]
            - int for polyalanine
        """
        if type(arg) is list and len(arg) > 0 and arg[0].__class__.__name__ is 'Atom':
            self.atoms = arg
        elif hasattr(arg, 'atoms') and type(arg.atoms) is list \
                and len(arg.atoms) > 0 and arg.atoms[0].__class__.__name__ is 'Atom':
            self.atoms = arg.atoms
        elif type(arg) is str:
            self.atoms = []
            if ':' in arg:
                seq, sec = arg.split(':')
                if len(sec) != len(seq):
                    raise Exception('Sequence length != secondary structure in ' + arg + ' !!!')
            else:
                seq = arg
                sec = 'C' * len(seq)

            for i, ch in enumerate(seq):
                self.atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i+1,
                        name=' CA ',
                        resname=aa_to_long(ch),
                        resnum=i+1,
                        occ=CABS_SS.get(sec[i], 1)
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
            self.atoms = []

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, index):
        return self.atoms[index]

    def __setitem__(self, index, atom):
        self.atoms[index] = atom

    def __delitem__(self, index):
        del self.atoms[index]

    def append(self, atom):
        self.atoms.append(atom)

    def extend(self, other):
        self.atoms.extend(other.atoms)

    def __repr__(self):
        return '\n'.join(str(atom) for atom in self.atoms)

    def __eq__(self, other):
        if len(self.atoms) != len(other.atoms):
            return False
        for atom1, atom2 in zip(self.atoms, other.atoms):
            if atom1 != atom2:
                return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def residues(self):
        """
        Returns a list of Atoms objects representing residues.
        """
        res = []
        residue = Atoms()
        residue.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_residue(residue.atoms[-1]):
                residue.append(atom)
            else:
                if residue not in res:
                    res.append(residue)
                for r in res:
                    if atom.same_residue(r[0]):
                        r.append(atom)
                        break
                else:
                    residue = Atoms()
                    residue.append(atom)
        if residue not in res:
            res.append(residue)
        return res

    def chains(self):
        """
        Returns a list of Atoms objects representing chains.
        """
        chn = []
        chain = Atoms()
        chain.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_chain(chain.atoms[-1]):
                chain.append(atom)
            else:
                if chain not in chn:
                    chn.append(chain)
                for ch in chn:
                    if atom.same_chain(ch[0]):
                        ch.append(atom)
                        break
                else:
                    chain = Atoms()
                    chain.append(atom)
        if chain not in chn:
            chn.append(chain)
        return chn

    def models(self):
        """
        Returns a list of Atoms objects representing models.
        """
        mdl = []
        model = Atoms()
        model.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_model(model.atoms[-1]):
                model.append(atom)
            else:
                if model not in mdl:
                    mdl.append(model)
                for m in mdl:
                    if atom.same_model(m[0]):
                        m.append(atom)
                        break
                else:
                    model = Atoms()
                    model.append(atom)
        if model not in mdl:
            mdl.append(model)
        return mdl

    def to_matrix(self):
        """"
        Returns numpy.matrix(N, 3) where N is number of Atoms.
        Matrix holds Atoms' coordinates.
        """
        return np.concatenate([a.coord.to_matrix() for a in self.atoms])

    def from_matrix(self, matrix):
        """
        Sets Atoms' coordinates from numpy.matrix(3,N) or(N,3).
        """
        if matrix.shape == (3, len(self)):
            for index, atom in enumerate(self.atoms):
                atom.coord = Vector3d(
                    matrix[0, index],
                    matrix[1, index],
                    matrix[2, index]
                )
        elif matrix.shape == (len(self), 3):
            for index, atom in enumerate(self.atoms):
                atom.coord = Vector3d(
                    matrix[index, 0],
                    matrix[index, 1],
                    matrix[index, 2]
                )
        else:
            raise Exception('Invalid matrix shape: ' + matrix.shape)
        return self

    def move(self, v):
        """
        Move atoms by vector.
        """
        if v is not None:
            for a in self.atoms:
                a.coord += v
        return self

    def rotate(self, matrix):
        """
        Rotate atoms by rotation matrix.
        """
        if matrix.shape != (3, 3):
            raise Exception('Invalid matrix shape: ' + matrix.shape)
        self.from_matrix(matrix * self.to_matrix().T)
        return self

    def cent_of_mass(self):
        """
        Returns a vector of the geometrical center of atoms.
        """
        com = Vector3d()
        for atom in self.atoms:
            com += atom.coord
        return com / len(self)

    def center_at_origin(self):
        """
        Moves atoms so that their geometrical center is in [0, 0, 0].
        """
        self.move(-self.cent_of_mass())
        return self

    def move_to(self, v):
        """
        Moves atoms so that their geometrical center is in [vx, vy, vz]. 
        """
        self.move(v - self.cent_of_mass())
        return self

    def compute_rotation(self, other, concentric=False):
        """
        Computes the rotation matrix for best fit between two sets of Atoms.
        """
        if len(self) != len(other):
            raise Exception('Atom sets have different length: %i != %i' % (len(self), len(other)))
        t = other.to_matrix()
        q = self.to_matrix()
        return kabsch(t, q, concentric)

    def str_align(self, other):
        """
        Aligns structurally set of atoms to another set.
        """
        r = self.compute_rotation(other)
        self.center_at_origin().rotate(r).move(other.cent_of_mass())
        return self

    def rmsd(self, other):
        """
        Calculates rmsd between two sets of atoms.
        """
        if len(self) != len(other):
            raise Exception('Atom sets have different length: %i != %i' % (len(self), len(other)))
        r = 0
        for a1, a2 in zip(self, other):
            r += a1.dist2(a2)
        return sqrt(r/len(self))

    def min_distance(self, other):
        """
        Calculates minimal distance between two sets of atoms.
        """
        return min(a.min_distance(other) for a in self.atoms)

    def change_chid(self, old, new):
        """
        Changes chain ID. Can do multiple changes at once.
        """
        if len(old) == len(new) and len(old) > 0:
            d = {}
            for o, n in zip(old, new):
                d[o] = n
            for atom in self:
                if atom.chid in d:
                    atom.chid = d[atom.chid]
        return self

    def model_count(self):
        """
        Returns number of models in Atoms object.
        """
        return len(self.models())

    def chain_count(self):
        """
        Returns number of chains in Atoms object.
        """
        return len(self.chains())

    def residue_count(self):
        """
        Returns number of residues in Atoms object.
        """
        return len(self.residues())

    def list_chains(self):
        """
        Returns a dictionary [chain ID] = chain_residue_count
        """
        d = {}
        for ch in self.chains():
            d[ch[0].chid] = Atoms(ch).residue_count()
        return d

    def max_dimension(self):
        """
        Returns maximal distance between any two atoms from the Atoms object.
        """
        return sqrt(max([p[0].dist2(p[1]) for p in combinations(self.atoms, 2)]))

    def make_pdb(self, bar_msg=''):
        """
        Returns a pdb-like formatted string.
        """
        models = self.models()
        if bar_msg:
            bar = ProgressBar(len(models), bar_msg)
        else:
            bar = None
        if len(models) == 1:
            s = self.__repr__()
        else:
            s = ''
            for m in models:
                s += 'MODEL%9i\n' % m[0].model
                s += m.__repr__()
                s += '\nENDMDL\n'
                if bar:
                    bar.update()
        if bar:
            bar.done(False)
        return s

    def save_to_pdb(self, filename, bar_msg=''):
        with open(filename, 'w') as f:
            f.write(self.make_pdb(bar_msg=bar_msg))

    def select(self, sele):
        """
        Selects subset of atoms defined by selection sentence.
        """
        if type(sele) is str:
            s = Selection(sele)
        else:
            s = sele
        return Atoms([a for a in self.atoms if a.match(s)])

    def drop(self, sele):
        """
        Removes subset of atoms defined by selection sentence.
        """
        if type(sele) is str:
            s = Selection(sele)
        else:
            s = sele
        return self.select(~s)

    def update_sec(self, sec):
        """
        Reads secondary structure dictionary sec[] with Atoms.resid_id() as keys
        and puts it into self.occ in CABS code:
        Helix - > 2.0, Sheet -> 4.0, Turn -> 3.0, Coil -> 1.0
        """
        if sec:
            for i in self.atoms:
                k = i.resid_id()
                if k in sec:
                    if sec[k] == 'H':
                        i.occ = 2.0
                    elif sec[k] == 'E':
                        i.occ = 4.0
                    elif sec[k] == 'T':
                        i.occ = 3.0
                    else:
                        i.occ = 1.0
        return self

    def update_bfac(self, bfac, d=0):
        for i in self.atoms:
            k = i.resid_id()
            if k in bfac:
                i.bfac = bfac[k]
            else:
                i.bfac = d

    def set_bfac(self, b=0.0):
        for atom in self.atoms:
            atom.bfac = b
        return self

    def valid_residues(self, must_have='CA, N, C, O'):
        """
        Returns only those residues that have atoms specified in "must_have" parameter.
        TODO: This is just temporary and it will be replaced by conditional selection class.
        """
        valid = Atoms()
        mh = [word.strip() for word in must_have.split(',')]
        for residue in self.residues():
            keep = True
            for nm in mh:
                if len(residue.select('name ' + nm)) == 0:
                    keep = False
                    break
            if keep:
                valid += residue
        return valid

    def remove_alternative_locations(self):
        """
        Removes atoms with alternative locations other than ' ' or 'A'
        """
        self.atoms = [atom for atom in self.atoms if (atom.alt == ' ' or atom.alt == 'A')]
        return self

    def set_model_number(self, number):
        """
        Sets model number to [number] for all atoms.
        """
        for a in self.atoms:
            a.model = number
        return self

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


class Selection:
    """
    Class representing atom selection.
    TODO: aliases
    """

    ARGS_KEYWORDS = ['MODEL', 'CHAIN', 'RESNUM', 'RESNAME', 'NAME']
    NO_ARGS_KEYWORDS = ['HETERO']
    OPERATORS = ['NOT', 'AND', 'OR']
    PARENTHESIS = ['(', ')']
    JOINTS = OPERATORS + PARENTHESIS
    KEYWORDS = ARGS_KEYWORDS + NO_ARGS_KEYWORDS + JOINTS

    def __init__(self, s=''):
        """
        Takes a string as input and parses it into selection tokens
        """
        self.tokens = []
        if s is not None:
            for word in s.replace('(', ' ( ').replace(')', ' ) ').split():
                if word.upper() in self.KEYWORDS:
                    self.tokens.append(word)
                elif len(self.tokens) is not 0:
                    self.tokens[-1] += " " + word

    def __invert__(self):
        return Selection("not( " + repr(self) + " )")

    def __repr__(self):
        return " ".join(self.tokens)


if __name__ == '__main__':
    pass

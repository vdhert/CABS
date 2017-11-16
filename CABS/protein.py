"""
Classes Protein, Peptide, ProteinComplex - prepares initial complex.
"""

import re
from copy import deepcopy
from random import randint
from string import ascii_uppercase
from math import exp

from CABS import utils
from CABS import logger
from CABS.pdblib import Pdb
from CABS.atom import Atoms
from CABS.vector3d import Vector3d

_name = 'Protein'


class Protein(Atoms):
    """
    Class for the protein molecule.
    """

    def __init__(self, source, flexibility=None, exclude=None, weights=None, work_dir='.'):

        Atoms.__init__(self)

        pdb = Pdb(source=source, selection='name CA')
        self.atoms = pdb.atoms.models()[0]
        logger.info(module_name=_name, msg="Loading %s as input protein" % source)

        # setup flexibility
        if flexibility:
            try:
                bfac = float(flexibility)
                self.atoms.set_bfac(bfac)
            except ValueError:
                if flexibility.lower() == 'bf':
                    pass
                elif flexibility.lower() == 'bfi':
                    for a in self.atoms:
                        if a.bfac > 1.:
                            a.bfac = 0.
                        elif a.bfac < 0.:
                            a.bfac = 1.
                        else:
                            a.bfac = 1. - a.bfac
                elif flexibility.lower() == 'bfg':
                    for a in self.atoms:
                        if a.bfac < 0.:
                            a.bfac = 1.
                        else:
                            a.bfac = exp(-0.5 * a.bfac ** 2)
                else:
                    try:
                        d, de = self.read_flexibility(flexibility)
                        self.atoms.update_bfac(d, de)
                    except IOError:
                        logger.warning(_name, 'Could not read flexibility file: %s' % flexibility)
                        logger.warning(_name, 'Using default flexibility(1.0) for all residues.')
                        self.atoms.set_bfac(1.0)
        else:
            self.atoms.set_bfac(1.0)

        # setup excluding
        self.exclude = {}
        if exclude:
            for s in exclude:
                words = s.split('@')
                if len(words) == 1:
                    key = 'ALL'
                else:
                    key = utils.pep2pep1(words[-1])
                if key in self.exclude:
                    self.exclude[key] += '+' + words[0]
                else:
                    self.exclude[key] = words[0]

            for k, v in self.exclude.items():
                self.exclude[k] = []
                for word in v.split('+'):
                    if ':' in word:
                        if '-' in word:
                            beg, end = word.split('-')
                            self.exclude[k].extend(self.atoms.atom_range(beg, end))
                        else:
                            self.exclude[k].append(word)
                    else:
                        chains = re.sub(r'[^%s]*' % word, '', ascii_uppercase)
                        self.exclude[k].extend(a.resid_id() for a in self.atoms.select('chain %s' % chains))

        ss = pdb.dssp(output=work_dir)
        self.old_ids = self.atoms.update_sec(ss).fix_broken_chains()
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for key, val in self.exclude.items():
            self.exclude[key] = [self.new_ids[r] for r in val]

        # setup rmsd_weights
        self.weights = None
        if weights and weights.lower() == 'flex':
            self.weights = [a.bfac for a in self.atoms]
        if weights and weights.lower() == 'ss':
            self.weights = [(a.occ + 1.) % 2 for a in self.atoms]
        else:
            try:
                default = 1.0
                self.weights = []
                weights_dict = {}
                with open(weights, 'rb') as _file:
                    for line in _file:
                        k, v = line.split()[:2]
                        weights_dict[k] = v

                if 'default' in weights_dict:
                    default = float(weights_dict['default'])

                for a in self.atoms:
                    w = weights_dict.get(a.resid_id())
                    w = float(w) if w else default
                    self.weights.append(w)

            except (IOError, ValueError):
                logger.warning(_name, 'Could not read weights file: %s' % weights)
                logger.warning(_name, 'Using default weights(1.0) for all atoms.')

        self.center = self.cent_of_mass()
        self.dimension = self.max_dimension()
        self.patches = {}

    def convert_patch(self, location):
        if location not in self.patches:
            chains = {}
            for res in [self.new_ids[r] for r in location.split('+')]:
                num, chid = res.split(':')
                if chid in chains:
                    chains[chid].append(num)
                else:
                    chains[chid] = [num]
            s = " or ".join(["(chain " + ch + " and resnum " + ",".join(chains[ch]) + ")" for ch in chains])
            patch = self.select(s)
            self.patches[location] = (patch.cent_of_mass() - self.center).norm()
        return self.patches[location]

    @staticmethod
    def read_flexibility(filename):

        key = r'[0-9A-Z]+:[A-Z]'
        val = r'[0-9.]+'

        patt_range = re.compile('(%s) *-* *(%s) +(%s)' % (key, key, val))
        patt_single = re.compile('(%s) +(%s)' % (key, val))

        with open(filename) as f:
            d = {}
            def_val = 1.0
            for line in f:
                if re.search('default', line):
                    def_val = float(line.split()[-1])
                else:
                    match = re.search(patt_range, line)
                    if match:
                        n1, c1 = match.group(1).split(':')
                        n2, c2 = match.group(2).split(':')
                        n1 = int(n1)
                        n2 = int(n2)
                        if c1 != c2 or n1 > n2:
                            raise Exception('Invalid range: \'%s\' in file: %s!!!' % (line, filename))
                        for i in range(n1, n2 + 1):
                            d[str(i) + ':' + c1] = float(match.group(3))
                    else:
                        match = re.search(patt_single, line)
                        if match:
                            d[match.group(1)] = float(match.group(2))
                        else:
                            raise Exception('Invalid syntax in flexibility file!!!')
            return d, def_val

    def generate_restraints(self, mode, gap, min_d, max_d):
        gap = int(gap)
        min_d = float(min_d)
        max_d = float(max_d)
        restr = []
        _len = len(self.atoms)

        for i in range(_len):
            a1 = self.atoms[i]
            ssi = int(a1.occ) % 2
            if mode == 'ss2' and ssi:
                continue
            for j in range(i + gap, _len):
                a2 = self.atoms[j]
                ssj = int(a2.occ) % 2
                if (mode == 'ss2' and ssj) or (mode == 'ss1' and ssi * ssj):
                    continue
                d = (a1.coord - a2.coord).length()
                if min_d < d < max_d:
                    if a1.bfac < a2.bfac:
                        w = a1.bfac
                    else:
                        w = a2.bfac
                    if w:
                        restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))
        return restr


class Peptide(Atoms):
    """
    Class for the peptides.
    """

    def __init__(self, source, conformation, location, work_dir='.'):
        logger.info(
            module_name=_name,
            msg='Loading ligand: {}, conformation - {}, location - {}'.format(
                source, conformation, location
            )
        )
        try:
            pdb = Pdb(source=source, selection='name CA', no_exit=True)
            atoms = pdb.atoms.models()[0]
            atoms.update_sec(pdb.dssp(output=work_dir))
        except Pdb.InvalidPdbInput:
            atoms = Atoms(source)
        atoms.set_bfac(0.0)
        self.conformation = conformation
        self.location = location
        Atoms.__init__(self, atoms)

    def random_conformation(self, lib=utils.RANDOM_LIGAND_LIBRARY):
        length = len(self)
        models, max_length, dim = lib.shape
        if length > max_length:
            raise Exception('Cannot generate random coordinates for peptide length = %d (max is %d)'
                            % (length, max_length))
        model = randint(0, models - 1)
        index = randint(0, max_length - length)
        self.from_numpy(lib[model][index: index + length])
        return self


class ProteinComplex(Atoms):
    """
    Class that assembles the initial complex.
    """

    def __init__(self, protein, flexibility, exclude, weights, peptides, replicas,
                 separation, insertion_attempts, insertion_clash, work_dir):
        logger.debug(module_name=_name, msg="Preparing the complex")
        Atoms.__init__(self)

        self.protein = Protein(
            protein, flexibility=flexibility, exclude=exclude, weights=weights, work_dir=work_dir
        )
        self.chain_list = self.protein.list_chains()
        self.protein_chains = ''.join(self.chain_list.keys())
        self.old_ids = deepcopy(self.protein.old_ids)

        self.peptides = []
        self.peptide_chains = ''
        if peptides:
            taken_chains = self.protein_chains + 'X'
            for num, p in enumerate(peptides):
                peptide = Peptide(*p, work_dir=work_dir)
                if peptide[0].chid in taken_chains:
                    peptide.change_chid(peptide[0].chid, utils.next_letter(taken_chains))
                taken_chains += peptide[0].chid
                self.peptide_chains += peptide[0].chid
                self.peptides.append(peptide)
                self.old_ids.update({atom.resid_id(): '%i:PEP%i' % (i + 1, num + 1) for i, atom in enumerate(peptide)})
                self.chain_list.update(peptide.list_chains())
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        exclude = []
        for key, value in self.protein.exclude.items():
            if key == 'ALL':
                kword = 'PEP'
            else:
                kword = key
            keys = [v for k, v in self.new_ids.items() if re.search(kword, k)]
            exclude.extend((r1, r2) for r1 in keys for r2 in value)
        self.protein.exclude = list(set(exclude))

        for i in range(replicas):
            model = deepcopy(self.protein)
            model.set_model_number(i + 1)
            for peptide in self.peptides:
                for attempt in range(insertion_attempts):
                    self.insert_peptide(self.protein, peptide, separation)
                    if model.min_distance(peptide) > insertion_clash:
                        peptide = deepcopy(peptide)
                        peptide.set_model_number(i + 1)
                        model.atoms.extend(peptide)
                        break
                else:
                    raise Exception('Maximum number of attempts to insert peptide %s reached!!!' % peptide.name)
            self.atoms.extend(model)
        logger.debug(module_name=_name, msg="Complex successfully created")

    @staticmethod
    def insert_peptide(protein, peptide, separation):

        radius = 0.5 * protein.dimension + separation

        if peptide.location == 'keep':
            location = peptide.cent_of_mass()
        elif peptide.location == 'random':
            location = Vector3d().random() * radius + protein.center
        else:
            location = protein.convert_patch(peptide.location) * radius + protein.center

        if peptide.conformation == 'random':
            peptide.random_conformation()

        peptide.move_to(location)


if __name__ == '__main__':
    pass

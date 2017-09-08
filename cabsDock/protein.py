"""
Classes Receptor, Ligand, Protein - prepares initial complex.
"""

import re
import logger
from copy import deepcopy
from os.path import exists, join
from random import randint
from string import ascii_uppercase

from cabsDock.atom import Atoms
from cabsDock.pdb import Pdb, InvalidPdbCode
from cabsDock.vector3d import Vector3d
from cabsDock.utils import AA_NAMES, RANDOM_LIGAND_LIBRARY, next_letter, fix_residue, check_peptide_sequence
from cabsDock.utils import PEPtoPEP1 as PP

__all__ = ["Protein"]
class Receptor(Atoms):
    """
    Class for the protein receptor molecule. Initialized with job's config dictionary.
    """

    def __init__(self, config):

        Atoms.__init__(self)

        name = config['receptor']
        selection = 'name CA and not HETERO'
        pdb = Pdb(name, selection=selection)
        self.atoms = pdb.atoms.models()[0]
        logger.info(module_name=__all__[0], msg = "Loading %s as receptor" % name)
        token = config.get('receptor_flexibility')
        if token:
            try:
                bfac = float(token)
                self.atoms.set_bfac(bfac)
            except ValueError:
                if token.lower() == 'bf':
                    pass
                elif token.lower() == 'bfi':
                    for a in self.atoms:
                        if a.bfac > 1.:
                            a.bfac = 0.
                        else:
                            a.bfac = 1. - a.bfac
                elif exists(token):
                    d, de = self.read_flexibility(token)
                    self.atoms.update_bfac(d, de)
                elif exists(join(config['work_dir'], token)):
                    d, de = self.read_flexibility(join(config['work_dir'], token))
                    self.atoms.update_bfac(d, de)
                else:
                    raise Exception('Invalid receptor_flexibility setting in \'%s\'!!!' % token)
        else:
            self.atoms.set_bfac(1.0)

        self.exclude = {}
        token = config.get('exclude')
        if token:
            for s in token:
                words = s.split('@')
                if len(words) == 1:
                    key = 'ALL'
                else:
                    key = PP(words[-1])
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

        self.old_ids = self.atoms.update_sec(pdb.dssp(dssp_command=config['dssp_command'],output=config['work_dir'])).fix_broken_chains()
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for key, val in self.exclude.items():
            self.exclude[key] = [self.new_ids[r] for r in val]

        self.center = self.cent_of_mass()
        self.dimension = self.max_dimension()
        self.patches = {}
        self.check_residue_modifications()

    def check_residue_modifications(self):
        for atom in self:
            atom.resname = fix_residue(atom.resname)
        return self

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
        l = len(self.atoms)

        for i in range(l):
            a1 = self.atoms[i]
            ssi = int(a1.occ) % 2
            if mode == 'ss2' and ssi:
                continue
            for j in range(i + gap + 1, l):
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


class Ligand(Atoms):
    """
    Class for the peptides.
    """

    def __init__(self, config, num):
        self.name, self.conformation, self.location = config['ligand'][num]
        self.selection = 'name CA and not HETERO'
        logger.info(module_name=__all__[0],
                    msg = "Loading ligand: name = %s, conformation = %s, location = %s" %
                          (self.name.split(':')[0], self.conformation, self.location) )
        try:
            pdb = Pdb(self.name, selection=self.selection)
            atoms = pdb.atoms.models()[0]
            atoms.update_sec(pdb.dssp(output=config['work_dir']))
        except InvalidPdbCode:
            logger.debug(module_name=__all__[0],msg = 'Provided ligand is not a valid pdb code/file')
            seq = self.name.split(':')[0]
            check_peptide_sequence(seq)
            atoms = Atoms(self.name)
        atoms.set_bfac(0.0)
        Atoms.__init__(self, atoms)
        # checks the input peptide sequence for non-standard amino acids.
        rev_dct = dict(map(reversed, AA_NAMES.items()))
        [check_peptide_sequence(rev_dct[peptide.resname]) for peptide in atoms.atoms]

    def random_conformation(self, lib=RANDOM_LIGAND_LIBRARY):
        length = len(self)
        models, max_length, dim = lib.shape
        if length > max_length:
            raise Exception('Cannot generate random coordinates for peptide length = %d (max is %d)'
                            % (length, max_length))
        model = randint(0, models - 1)
        index = randint(0, max_length - length)
        self.from_matrix(lib[model][index: index + length])
        return self


class ProteinComplex(Atoms):
    """
    Class that assembles the initial complex.
    """

    def __init__(self, config):
        logger.debug(module_name=__all__[0], msg = "Preparing the complex")
        Atoms.__init__(self)
        self.separation = config['initial_separation']

        receptor = Receptor(config)
        self.chain_list = receptor.list_chains()
        self.receptor_chains = ''.join(self.chain_list.keys())
        self.old_ids = deepcopy(receptor.old_ids)

        ligands = []
        self.ligand_chains = ''
        if 'ligand' in config:
            taken_chains = self.receptor_chains + 'X'
            for num, ligand in enumerate(config['ligand']):
                l = Ligand(config, num)
                if l[0].chid in taken_chains:
                    l.change_chid(l[0].chid, next_letter(taken_chains))
                taken_chains += l[0].chid
                self.ligand_chains += l[0].chid
                ligands.append(l)
                self.old_ids.update({atom.resid_id(): '%i:PEP%i' % (i + 1, num + 1) for i, atom in enumerate(l)})
                self.chain_list.update(l.list_chains())
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        exclude = []
        for key, value in receptor.exclude.items():
            if key == 'ALL':
                kword = 'PEP'
            else:
                kword = key
            keys = [v for k, v in self.new_ids.items() if re.search(kword, k)]
            exclude.extend((r1, r2) for r1 in keys for r2 in value)
        receptor.exclude = list(set(exclude))

        for i in range(config['replicas']):
            model = deepcopy(receptor)
            model.set_model_number(i + 1)
            for ligand in ligands:
                for attempt in range(config['ligand_insertion_attempts']):
                    self.insert_ligand(receptor, ligand)
                    if model.min_distance(ligand) > config['ligand_insertion_clash']:
                        ligand = deepcopy(ligand)
                        ligand.set_model_number(i + 1)
                        model.atoms.extend(ligand)
                        break
                else:
                    raise Exception('Maximum number of attempts to insert ligand %s reached!!!' % ligand.name)
            self.atoms.extend(model)
            self.receptor = receptor
            self.ligands = ligands
        logger.debug(module_name=__all__[0], msg="Complex successfully created")

    def insert_ligand(self, receptor, ligand):

        radius = 0.5 * receptor.dimension + self.separation

        if ligand.location == 'keep':
            location = ligand.cent_of_mass()
        elif ligand.location == 'random':
            location = Vector3d().random() * radius + receptor.center
        else:
            location = receptor.convert_patch(ligand.location) * radius + receptor.center

        if ligand.conformation == 'random':
            ligand.random_conformation()

        ligand.move_to(location)

if __name__ == '__main__':
    pass

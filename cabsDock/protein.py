"""
Classes Receptor, Ligand, Protein - prepares initial complex.
"""

import re
from copy import deepcopy
from os.path import exists, join, isfile
from random import randint

from cabsDock.utils import check_peptide_sequence
from cabsDock.utils import AA_NAMES

from atom import Atoms
from pdb import Pdb, InvalidPdbCode
from utils import RANDOM_LIGAND_LIBRARY, next_letter, fix_residue
from vector3d import Vector3d


class Receptor(Atoms):
    """
    Class for the protein receptor molecule. Initialized with job's config dictionary.
    """

    def __init__(self, config):
        name = config['receptor']
        selection = 'name CA and not HETERO'
        if isfile(name):
            pdb = Pdb(pdb_file=name)
        elif isfile(join(config['work_dir'], name)):
            pdb = Pdb(pdb_file=join(config['work_dir'], name))
        else:
            pdb = Pdb(pdb_code=name[:4])
            m = re.match(r'.{4}:([A-Z]*)', name)
            if m:
                selection += ' and chain ' + ','.join(m.group(1))
        atoms = pdb.atoms.remove_alternative_locations().select(selection).models()[0]

        if 'receptor_flexibility' in config:
            token = config['receptor_flexibility']
            try:
                bfac = float(token)
                atoms.set_bfac(bfac)
            except ValueError:
                if token.lower() == 'bfac':
                    pass
                elif token.lower() == 'bfac_inv':
                    for a in atoms:
                        a.bfac = 1. - 0.01 * a.bfac
                elif exists(token):
                    d, de = self.read_flexibility(token)
                    atoms.update_bfac(d, de)
                elif exists(join(config['work_dir'], token)):
                    d, de = self.read_flexibility(join(config['work_dir'], token))
                    atoms.update_bfac(d, de)
                else:
                    raise Exception('Invalid receptor_flexibility setting in \'%s\'!!!' % token)
        else:
            atoms.set_bfac(1)

        self.old_ids = atoms.update_sec(pdb.dssp(dssp_command=config['dssp_command'])).fix_broken_chains()
        # self.new_ids = {v: k for k, v in self.old_ids.items()}
        Atoms.__init__(self, atoms)
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
        patt = re.compile(r'(.*)WEIGHT(.*)', re.IGNORECASE)

        with open(filename) as f:
            d = {}
            def_val = None
            for line in f:
                if not def_val:
                    def_val = float(line)
                match = re.match(patt, line)
                if match:
                    key = match.group(1).strip()
                    val = float(match.group(2))
                    k1, k2 = key.split('-')
                    n1, c1 = k1.strip().split(':')
                    n2, c2 = k2.strip().split(':')
                    n1 = int(n1)
                    n2 = int(n2)
                    if c1 != c2 or n1 > n2:
                        raise Exception('Invalid range: \'%s\' in file: %s!!!' % (key, filename))
                    for i in range(n1, n2 + 1):
                        d[str(i) + ':' + c1] = val
            return d, def_val

    def generate_restraints(self, gap, min_d, max_d):
        restr = []
        l = len(self.atoms)

        for i in range(l):
            a1 = self.atoms[i]
            for j in range(i + gap + 1, l):
                a2 = self.atoms[j]
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
        selection = 'name CA and not HETERO'
        if exists(self.name):
            pdb = Pdb(pdb_file=self.name)
            atoms = pdb.atoms.remove_alternative_locations().select(selection).models()[0]
            atoms.update_sec(pdb.dssp())
        elif exists(join(config['work_dir'], self.name)):
            pdb = Pdb(pdb_file=join(config['work_dir'], self.name))
            atoms = pdb.atoms.remove_alternative_locations().select(selection).models()[0]
            atoms.update_sec(pdb.dssp())
        else:
            try:
                pdb = Pdb(pdb_code=self.name[:4])
                m = re.match(r'.{4}:([A-Z]*)', self.name)
                if m:
                    selection += ' and chain ' + ','.join(m.group(1))
                atoms = pdb.atoms.remove_alternative_locations().select(selection).models()[0]
                atoms.update_sec(pdb.dssp())
            except InvalidPdbCode:
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
                print(ligand)
                ligands.append(l)
                self.old_ids.update({atom.resid_id(): '%i:LIG%i' % (i + 1, num + 1) for i, atom in enumerate(l)})
                self.chain_list.update(l.list_chains())
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for i in range(config['replicas']):
            model = deepcopy(receptor)
            model.set_model_number(i + 1)
            for ligand in ligands:
                for attempt in range(config['ligand_insertion_attempts']):
                    self.insert_ligand(receptor, ligand)
                    if model.min_distance(ligand) > config['ligand_insertion_clash']:
                        ligand = deepcopy(ligand)
                        ligand.set_model_number(i + 1)
                        model.extend(ligand)
                        break
                else:
                    raise Exception('Maximum number of attempts to insert ligand %s reached!!!' % ligand.name)
            self.atoms.extend(model)
            self.receptor = receptor
            self.ligands = ligands

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

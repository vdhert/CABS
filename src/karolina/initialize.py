from abc import ABCMeta, abstractmethod

class Initialize(object):

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self):
        pass

    def initialize_residues(self):
        from residue import Residue

        if isinstance(self, Residue):
            raise TypeError('Cannot call initialize_residues() method on a Residue object.')
        else:
            atoms_residues = (self.atoms).residues()
            residues = list()
            for atoms in atoms_residues:
                residues.append(Residue(atoms, (atoms.atoms[0]).resname, (atoms.atoms[0]).resnum))
            return residues

    def initialize_chains(self):
        from residue import Residue
        from chain import Chain

        if isinstance(self, Residue) or isinstance(self, Chain):
            raise TypeError('Cannot call initialize_chains() method on either Residue or Chain object.')
        else:
            atoms_chains = (self.atoms).chains()
            chains = list()
            for atoms in atoms_chains:
                chains.append(Chain(atoms, (atoms.atoms[0]).chid))
            return chains

    def initialize_models(self, headers, num):
        from residue import Residue
        from chain import Chain
        from model import Model

        if isinstance(self, Residue) or isinstance(self, Chain) or isinstance(self, Model):
            raise TypeError('Cannot call initialize_models() method on Residue, Chain or Model object.')
        else:
            atoms_models = (self.atoms).models()
            models = list()
            for atoms in atoms_models:
                models.append(Model(atoms, headers[atoms_models.index(atoms)], (atoms.atoms[0]).model, num))
            return models

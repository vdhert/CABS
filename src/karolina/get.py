from abc import ABCMeta, abstractmethod
from initialize import Initialize

class Get(Initialize):

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self):
        pass

    def get_atoms(self):
        return self.atoms

    def get_residues(self):
        from residue import Residue

        if isinstance(self, Residue):
            raise TypeError('Cannot call get_residues() method on a Residue object.')
        elif hasattr(self, 'residues'):
            return self.residues
        else:
            return self.initialize_residues()

    def get_chains(self, **kwargs):

        '''
        If 'chain' arg is provided, returns only chains with chain ids specified in the argument.
        '''

        from trajectory import Trajectory
        from chain import Chain
        from residue import Residue

        if isinstance(self, Chain) or isinstance(self, Residue):
            raise TypeError('Cannot call get_chains() method on either Chain or Residue object.')
        if kwargs:
            if 'chain' in kwargs:
                chains = list()
                if isinstance(self, Trajectory):
                    for chid in kwargs['chain']:
                        for model in self.models:
                            chains.append(Chain((model.atoms).select('CHAIN ' + str(chid)), chid))
                else:
                    for chid in kwargs['chain']:
                        chains.append(Chain((self.atoms).select('CHAIN ' + str(chid)), chid))
                return chains
            else:
                raise KeyError('Wrong name of the argument.')
        else:
            if hasattr(self, 'chains'):
                return self.chains
            else:
                return self.initialize_chains()

    def get_models(self):
        from residue import Residue
        from chain import Chain
        from model import Model

        if isinstance(self, Residue) or isinstance(self, Chain) or isinstance(self, Model):
            raise TypeError('Cannot call get_models() method on Residue, Chain or Model object.')
        elif hasattr(self, 'models'):
            return self.models
        else:
            return self.initialize_models()

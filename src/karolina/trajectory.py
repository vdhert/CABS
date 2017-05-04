from atom import Atoms
from model import Model
from get import Get

class Trajectory(Get):

    def __init__(self, atoms, receptor_id, ligand_id, num, headers, **kwargs):

        '''
        -- 'atoms' arg is an Atoms object
        -- 'receptor_id' is a string
        -- 'ligand_id' is a list containing strings
        -- 'num' is an integer
        -- 'headers' is a list containing lists, which contain floats
        '''

        self.models = list()
        if isinstance(atoms, Atoms):
            self.atoms = atoms
            self.models = self.initialize_models(headers, num)
        elif isinstance(atoms, list):
            for i, j in zip(atoms, range(1, 1001)):
                if isinstance(i, Model):
                    (i.atoms).set_model_number(j)
                    i.num = j
                    (self.models).append(i)
                else:
                    raise TypeError('Cannot initialize a Trajectory object from a list containing objects other than a Model object.')
        else:
            raise TypeError('Initialization of a Trajectory object only possible from an Atoms object or a list containing Model objects.')
        self.num = num
        self.receptor_id = receptor_id
        self.ligand_id = ligand_id
        self.__dict__.update()

    def __iter__(self):
        for model in self.models:
            yield model

    def __getitem__(self, index):
        for i in range(len(self.models)):
            if i == index:
                return self.models[i]

    def __len__(self):
        return len(self.models)

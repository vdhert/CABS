from get import Get

class Residue(Get):

    def __init__(self, atoms, resname, resnum, **kwargs):

        '''
        -- 'atoms' arg is an Atoms object
        -- 'resname' arg is a string
        -- 'resnum' arg is an integer
        '''

        self.atoms = atoms
        self.resname = resname
        self.resnum = resnum
        self.__dict__.update()

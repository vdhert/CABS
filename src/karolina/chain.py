from get import Get

class Chain(Get):

    def __init__(self, atoms, chid, **kwargs):

        '''
        -- 'atoms' arg is an Atoms object
        -- 'chid' arg is a string
        '''

        self.atoms = atoms
        self.chid = chid
        self.residues = self.initialize_residues()
        self.__dict__.update()

from get import Get

class Model(Get):

    def __init__(self, atoms, header, num, traj_num, **kwargs):

        '''
        -- 'atoms' arg is an Atoms object
        -- 'header' arg is a list containing floats
        -- 'num' arg is an integer
        -- 'traj_num' is an integer
        '''

        self.atoms = atoms
        self.num = num
        self.chains = self.initialize_chains()
        self.temp = header[0]
        self.e_rec = header[1]
        self.e_lig = header[2]
        self.e_int = header[3]
        self.e_tot = header[4]
        self.traj_num = traj_num
        self.__dict__.update()

    def __iter__(self):
        for chain in self.chains:
            yield chain

    def __getitem__(self, index):
        for i in range(len(self.chains)):
            if i == index:
                return self.chains[i]

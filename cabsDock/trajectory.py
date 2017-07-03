from copy import deepcopy

import StringIO
import numpy
import numpy as np

from atom import Atom, Atoms
from pdb import Pdb
from utils import ranges, kabsch, ProgressBar
import warnings

__all__ = ['Trajectory', 'Header']


class Header:
    """Trajectory header read from CABS output: energies and temperatures"""

    class CannotMerge(Exception):
        """Raised when trying to merge headers for different frames"""

        def __init__(self, h1, h2):
            self.msg = 'Cannot merge headers: %s and %s' % (h1, h2)

        def __str__(self):
            return self.msg

    def __init__(self, line):
        header = line.split()
        self.model = int(header[0])
        self.length = (int(header[1]) - 2,)
        self.energy = np.matrix(header[2: -2], float)
        self.temperature = float(header[-2])
        self.replica = int(header[-1])
        self.rmsd = 0

    def __repr__(self):
        return 'Replica: %d Model: %d Length: %s T: %.2f E: %s' % (
            self.replica,
            self.model,
            self.length,
            self.temperature,
            str(self.energy.tolist())
        )

    def __add__(self, other):
        """Merges two headers from two chains of the same frame"""
        if self.replica != other.replica or self.model != other.model:
            raise Header.CannotMerge(self, other)
        else:
            dt = self.temperature - other.temperature
            if dt ** 2 > 1e-6:
                raise Exception("Cannot merge headers with different T!!!")
            else:
                h = deepcopy(self)
                h.length += other.length
                h.energy = np.concatenate([self.energy, other.energy])
        return h

    def get_energy(self, number_of_peptides):
        #number_of_peptides fixes energy calculations
        """ Returns the total energy of the system. """
        if not number_of_peptides:
            print("Unknown number of peptides. Assuming 1.")
            num_pept = 1
        else:
            num_pept = number_of_peptides
        int_submtrx_size = self.energy.shape[0]-num_pept
        int_enrg = np.sum(self.energy[:int_submtrx_size,-num_pept:])
        # print('whole')
        # print(self.energy)
        # print('int')
        # print(int_submtrx_size)
        # print (self.energy[:int_submtrx_size,-num_pept:])
        # print int_enrg
        return int_enrg


class Trajectory(object):
    """
    Class holding compressed trajectory.
    """
    GRID = 0.61

    def __init__(self, template, coordinates, headers):
        self.template = template
        self.coordinates = coordinates
        self.headers = headers
        self.rmsd_native = None
        self.number_of_peptides = None

    @staticmethod
    def read_seq(filename):
        atoms = []
        with open(filename) as f:
            for i, line in enumerate(f):
                atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i + 1,
                        name='CA',
                        alt=line[7],
                        resname=line[8:11],
                        chid=line[12],
                        resnum=int(line[1:5]),
                        icode=line[5],
                        occ=float(line[15]),
                        bfac=float(line[16:22])
                    )
                )
        return atoms

    @staticmethod
    def read_traf(filename):
        headers = []
        replicas = {}

        def save_header(h):
            headers.append(h)

        def save_coord(c, r):
            if r not in replicas:
                replicas[r] = []
            replicas[r].extend(c[3:-3])

        with open(filename) as f:
            current_header = None
            current_coord = []
            for line in f:
                if '.' in line:
                    header = Header(line)
                    if not current_header:
                        current_header = header
                    else:
                        save_coord(current_coord, current_header.replica)
                        current_coord = []
                        try:
                            current_header = current_header + header
                        except Header.CannotMerge:
                            save_header(current_header)
                            current_header = header
                else:
                    current_coord.extend(map(int, line.split()))
            save_header(current_header)
            save_coord(current_coord, current_header.replica)

        headers.sort(key=lambda x: x.model)
        headers.sort(key=lambda x: x.replica)
        coordinates = np.array([Trajectory.GRID * x for y in sorted(replicas) for x in replicas[y]])

        return headers, coordinates

    @classmethod
    def read_trajectory(cls, traf, seq):
        template = Atoms(Trajectory.read_seq(seq))
        headers, coordinates = Trajectory.read_traf(traf)

        replicas = len(set(h.replica for h in headers))
        if len(headers) % replicas:
            raise Exception('Replicas have different sizes!!!')
        models = len(headers) / replicas
        length = headers[0].length

        if any(length != h.length for h in headers):  # check if all frames have the same shape
            raise Exception('Invalid headers in %s!!!' % traf)

        sum_length = sum(length)

        if sum_length != len(template):
            # check if number of atoms in SEQ matches that in trajectory headers
            raise Exception('Different number of atoms in %s and %s!!!' % (traf, seq))

        # final test if information from headers agrees with number of coordinates
        size_test = replicas * models * sum_length * 3
        if size_test != len(coordinates):
            raise Exception('Invalid number of atoms in %s!!!' % traf)

        coordinates = coordinates.reshape(replicas, models, sum_length, 3)
        return cls(template, coordinates, headers)

    def select(self, selection):
        template = self.template.select(selection)
        inds = [self.template.atoms.index(a) for a in template]
        return Trajectory(template, self.coordinates[:,:,inds,:], self.headers)

    def to_atoms(self):
        result = Atoms()
        num = 0
        shape = self.coordinates.shape
        for model in self.coordinates.reshape(-1, len(self.template), 3):
            atoms = deepcopy(self.template)
            num += 1
            atoms.set_model_number(num)
            atoms.from_matrix(model)
            result.extend(atoms)
        self.coordinates.reshape(shape)
        return result

    def align_to(self, target, selection=''):
        if not selection:
            selection = 'chain ' + ','.join(target.list_chains().keys())
        pieces = ranges([self.template.atoms.index(a) for a in self.template.select(selection)])

        t = target.to_matrix()
        t_com = np.average(t, 0)
        t = np.subtract(t, t_com)

        shape = self.coordinates.shape
        for model in self.coordinates.reshape(-1, len(self.template), 3):
            query = np.concatenate([model[piece[0]:piece[1]] for piece in pieces])
            q_com = np.average(query, 0)
            q = np.subtract(query, q_com)
            np.copyto(model, np.add(np.dot(np.subtract(model, q_com), kabsch(t, q, concentric=True)), t_com))
        self.coordinates.reshape(shape)

    def rmsd_matrix(self, msg=''):
        """
        Calculates rmsd matrix with no fitting for all pairs od models in trajectory.
        :return: np.array
        """

        def rmsd(m1, m2, ml):
            return np.sqrt(np.sum((m1 - m2) ** 2) / ml)

        model_length = len(self.template)
        models = self.coordinates.reshape(-1, model_length, 3)
        dim = len(models)
        result = np.zeros((dim, dim))
        if msg:
            bar = ProgressBar((dim * dim - dim) / 2, msg=msg)
        else:
            bar = None
        for i in range(dim):
            for j in range(i + 1, dim):
                if bar:
                    bar.update()
                result[i, j] = result[j, i] = rmsd(models[i], models[j], model_length)
        if bar:
            bar.done(True)
        return result

    def rmsd_to_native(self, native_pdb="", native_receptor_chain="", native_peptide_chain="", model_peptide_chain=""):
        """
        Calculates a list of ligand - rmsd of the models to the native structure (argument 'native' is either
        a PDB code (to be downloaded) or a local PDB file).
        :return: np.array
        """

        def rmsd(m1, m2, length):
            return np.sqrt(np.sum((m1 - m2) ** 2) / length)

        target_selection = 'name CA and not HETERO'
        target_selection += ' and chain ' + ','.join(native_receptor_chain)
        pdb = Pdb(pdb_code=native_pdb[:4])
        native = pdb.atoms.remove_alternative_locations().select(target_selection).models()[0]
        shape = self.coordinates.shape
        self.align_to(native, target_selection)
        models_peptide_traj = self.select("chain " + model_peptide_chain)
        peptide_length = len(models_peptide_traj.template)
        models_peptide = models_peptide_traj.coordinates.reshape(-1, peptide_length, 3)
        # TO BE FIXED: cabsDock.atoms.to_matrix() method returns numpy.matrix, for consistency
        # should return numpy.array instead.
        native_peptide = numpy.array(pdb.atoms.remove_alternative_locations().select(
            "name CA and not HETERO and chain " + native_peptide_chain
        ).models()[0].to_matrix())
        result = np.zeros(
            (len(models_peptide))
        )
        for i, h in zip(range(len(models_peptide)), self.headers):
            result[i] = rmsd(models_peptide[i], native_peptide, peptide_length)
            h.rmsd = result[i]
        self.coordinates.reshape(shape)
        print('... done.')
        return result

    def get_model(self, model):
        """
        Do poprawy 
        """
        shape = self.coordinates.shape
        coordinates = self.coordinates.reshape(-1, len(self.template), 3)[model]
        atoms = deepcopy(self.template)
        atoms.set_model_number(model + 1)
        m = atoms.from_matrix(coordinates)
        self.coordinates.reshape(shape)
        return m

    def to_pdb(self, name = None, mode='models', to_dir = None):
        """
        Method for transforming a trajectory instance into a PDB file-like object.
        :param mode:    'models' -- the method returns a list of StringIO objects, each representing one model from the trajectory;
                        'replicas' -- the method returns a list of StringIO objects, each representing one replica from the trajectory.
        :return: StringIO object
        """
        execution_mode = {'models': (self.coordinates[0], 'model'), 'replicas': (self.coordinates, 'replica')}
        if to_dir:
            for i, m in enumerate(execution_mode[mode][0]):
                Trajectory(self.template, m, None).to_atoms().save_to_pdb(
                    (execution_mode[mode][1] if name is None else name)
                    +
                    ('' if len(execution_mode[mode][0]) == 1 else '_{0}'.format(i))
                    +
                    '.pdb'
                    )
            out = True
        else:
            out =  [
                StringIO.StringIO(
                    Trajectory(self.template, m, None).to_atoms().make_pdb()
                    )
                for m in execution_mode[mode][0]
                ]
        return out

    def rmsf(self, chains = ''):
        mdls = self.select('chain ' + ','.join(chains))
        mdls.align_to(mdls.get_model(1), 'chain ' + ','.join(chains))
        mdl_lth = len(mdls.template)
        mdls_crds = numpy.stack(mdls.coordinates.reshape(-1, mdl_lth, 3), axis=1)
        avg = [ numpy.mean(rsd, axis=0) for rsd in mdls_crds ]
        return [numpy.mean([numpy.linalg.norm(avg[i] - case) for case in rsd]) for i, rsd in enumerate(mdls_crds)]

if __name__ == '__main__':
    tra = Trajectory.read_trajectory('/Users/maciek/Desktop/1dkx/TRAF', '/Users/maciek/Desktop/1dkx/SEQ')
    tra.coordinates = tra.coordinates[0:1]
    print(tra.coordinates.shape)
    tra.to_pdb(name='blagh.pdb', mode='replicas', to_dir='/Users/maciek/Desktop')
    rmsf=tra.rmsf(chains='A')
    print numpy.mean(rmsf)
    from matplotlib import pyplot
    pyplot.plot(rmsf)
    pyplot.show()
from copy import deepcopy

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
        self.energy = np.matrix(header[4:1:-1], float)
        self.temperature = float(header[5])
        self.replica = int(header[6])
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
            de = self.energy[0, 0] - other.energy[0, 0]
            if dt ** 2 > 1e-6 or de ** 2 > 1e-6:
                raise Exception("Cannot merge headers with different T or E!!!")
            else:
                h = deepcopy(self)
                h.length += other.length
                h.energy = np.concatenate([self.energy, other.energy])
        return h


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

    @staticmethod
    def read_trajectory(traf, seq):
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
            # check if number of atoms in SEQ matches that in trajectory head.ers
            raise Exception('Different number of atoms in %s and %s!!!' % (traf, seq))

        # final test if information from headers agrees with number of coordinates
        size_test = replicas * models * sum_length * 3
        if size_test != len(coordinates):
            raise Exception('Invalid number of atoms in %s!!!' % traf)

        coordinates = coordinates.reshape(replicas, models, sum_length, 3)
        return Trajectory(template, coordinates, headers)

    def select(self, selection):
        template = self.template.select(selection)
        pieces = ranges([self.template.atoms.index(a) for a in template])
        shape = self.coordinates.shape
        self.coordinates.reshape(-1, len(self.template), 3)
        coordinates = []
        for replica in self.coordinates:
            replica_new = []
            for model in replica:
                replica_new.append(np.concatenate(
                    [
                        model[piece[0]:piece[1]] for piece in pieces
                    ]
                ))
            coordinates.append(np.stack(replica_new))
        coordinates = np.stack(coordinates)
        self.coordinates.reshape(shape)
        if len(self.coordinates) == 0:
            warnings.warn(
                "The selection \"{0}\" results in an empty trajectory.".format(
                    selection), UserWarning
                )
        return Trajectory(template, coordinates, self.headers)

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

    def filter(self, number):
        """
        Temporary filtering - top N from each replica by total energy
        """
        replicas = len(self.coordinates)
        models = len(self.coordinates[0])
        if 0 < replicas * models <= number:
            return self
        else:
            headers = []
            coordinates = []
            from_replica = int(number / replicas)
            remains = number - replicas * from_replica
            for replica in range(replicas):
                current = [h for h in self.headers if h.replica == replica + 1]
                by_energy = sorted(current, key=lambda x: np.sum(x.energy[:, 1:2]))
                top = from_replica + (replica < remains)
                current = sorted(sorted(by_energy[:top], key=lambda x: x.model), key=lambda x: x.replica)
                headers.extend(current)
                mdls = []
                for h in current:
                    mdls.append(h.model - 1)
                rngs = ranges(mdls)
                for r in rngs:
                    c = self.coordinates[replica][r[0]: r[1]]
                    coordinates.append(c)
            return Trajectory(self.template, np.concatenate(coordinates), headers)

    def rmsd_matrix(self, msg=''):
        """
        Calculates rmsd matrix with no fitting for all pairs od models in trajectory.
        :return: np.array
        """
        model_length = len(self.template)

        def rmsd(m1, m2):
            return np.sqrt(np.sum((m1 - m2) ** 2) / model_length)

        shape = self.coordinates.shape
        models = self.coordinates.reshape(-1, model_length, 3)
        dim = len(models)
        result = np.zeros(dim * dim).reshape(dim, dim)
        if msg:
            bar = ProgressBar((dim * dim - dim) / 2, msg=msg)
        else:
            bar = None
        for i in range(dim):
            for j in range(i + 1, dim):
                if bar:
                    bar.update()
                result[i, j] = result[j, i] = rmsd(models[i], models[j])
        if bar:
            bar.done(True)
        self.coordinates.reshape(shape)
        return result

    def rmsd_to_native(self, native_pdb="", native_receptor_chain="", native_peptide_chain="", model_peptide_chain=""):
        """
        Calculates a list of ligand - rmsd of the models to the native structure (argument 'native' is either
        a PDB code (to be downloaded) or a local PDB file).
        :return: np.array
        """

        def rmsd(m1, m2, length):
            return np.sqrt(np.sum(m1 - m2)**2 / length)

        target_selection = 'name CA and not HETERO'
        target_selection += ' and chain ' + ','.join(native_receptor_chain)
        pdb = Pdb(pdb_code=native_pdb[:4])
        native = pdb.atoms.remove_alternative_locations().select(target_selection).models()[0]
        shape = self.coordinates.shape
        self.align_to(native, target_selection)
        models_peptide_traj = self.select("chain " + model_peptide_chain)
        peptide_length = len(models_peptide_traj.template)
        models_peptide = models_peptide_traj.coordinates.reshape(-1, peptide_length, 3)
        native_peptide = pdb.atoms.remove_alternative_locations().select(
            "name CA and not HETERO and chain " + native_peptide_chain
        ).models()[0].to_matrix()
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

if __name__ == '__main__':
    tra = Trajectory.read_trajectory('.CABS/TRAF', '.CABS/SEQ')
    tra.rmsd_matrix()
    #tra.rmsd_to_native(native_pdb='1jbu', native_receptor_chain='H', native_peptide_chain ='X', model_peptide_chain='C')
    # from pdb import Pdb
    # target = Pdb(pdb_code='1jbu').atoms.select('name CA and chain H')
    # tra.align_to(target, 'chain H')
#     traf = tra.select('chain B, D')
#     traf.to_atoms().save_to_pdb('dupa.pdb')

import numpy as np
from copy import deepcopy

from atom import Atom, Atoms
from utils import ranges, kabsch

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


class Trajectory:
    """
    Class holding compressed trajectory.
    """
    GRID = 0.61

    def __init__(self, template, coordinates, headers):
        self.template = template
        self.coordinates = coordinates
        self.headers = headers

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
            # check if number of atoms in SEQ matches that in trajectory head.ers
            raise Exception('Different number of atoms in %s and %s!!!' % (traf, seq))

        # final test if information from headers agrees with number of coordinates
        size_test = replicas * models * sum_length * 3
        if size_test != len(coordinates):
            raise Exception('Invalid number of atoms in %s!!!' % traf)

        coordinates = coordinates.reshape(replicas, models, sum_length, 3)
        return cls(template, coordinates, headers)

    def select(self, selection):
        template = self.template.select(selection)
        return ranges([self.template.atoms.index(a) for a in template])

    def to_atoms(self):
        replicas, models = self.coordinates.shape[:2]
        result = Atoms()
        num = 0
        for r in range(replicas):
            for m in range(models):
                atoms = deepcopy(self.template)
                tail = '%d_%d' % (r + 1, m + 1)
                for a in atoms:
                    a.tail = tail
                num += 1
                atoms.set_model_number(num)
                atoms.from_matrix(self.coordinates[r][m])
                result.extend(atoms)
        return result

    def align_to(self, target, selection=''):
        if not selection:
            selection = 'chain ' + ','.join(target.list_chains().keys())
        rng = self.select(selection)
        if not rng:
            raise Exception('Invalid selection: \'%s\'' % selection)
        else:
            t = target.to_matrix()
            t_com = np.average(t, 0)
            t = np.subtract(t, t_com)

            shape = self.coordinates.shape
            if len(shape) == 4:
                models = shape[0] * shape[1]
            elif len(shape) == 3:
                models = shape[0]
            else:
                raise Exception('Invalid trajectory !!!')

            for m in range(models):
                model = self.coordinates[]
                sele = np.concatenate([model[r[0]: r[1]] for r in rng])
                sele_com = np.average(sele, 0)
                q = np.subtract(sele, sele_com)
                np.copyto(
                    model, np.add(np.dot(np.subtract(model, sele_com), kabsch(t, q, concentric=True)), t_com)
                )

    def filter(self, number):
        """
        Temporary filtering - top N from each replica by total energy
        """
        replicas, models = self.coordinates.shape[0:2]
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

if __name__ == '__main__':
    tra = Trajectory.read_trajectory('CABS/TRAF', 'CABS/SEQ')
    tra.align_to()

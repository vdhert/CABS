import numpy as np
from copy import deepcopy

from atom import Atom, Atoms
from utils import ProgressBar, ranges, line_count

__all__ = ['Trajectory', 'Header']


class Trajectory:
    def __init__(self, template, coordinates):
        self.template = template
        self.coordinates = coordinates

    def select(self, selection):
        template = self.template.select(selection)
        rng = ranges([self.template.atoms.index(a) for a in template])
        return Trajectory(template, np.concatenate([self.coordinates[:, :, r[0]: r[1]] for r in rng], 2))

    def to_atoms(self):
        replicas = []
        for r in range(self.coordinates.shape[0]):
            replicas.append(Atoms(None))
            for m in range(self.coordinates.shape[1]):
                atoms = deepcopy(self.template)
                atoms.set_model_number(m + 1)
                atoms.from_matrix(0.61* self.coordinates[r][m])
                replicas[r].extend(atoms)
        return replicas


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
        self.e_complex = None

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


def read_trajectory(traf, seq):
    template = Atoms(read_seq(seq))
    headers, coordinates = read_traf(traf)

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

    shape = (replicas, models, sum_length)
    coordinates = coordinates.reshape(shape + (3,))
    return (template, coordinates), headers


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


def read_traf(filename):
    headers = []
    replicas = {}

    def save_header(h):
        headers.append(h)

    def save_coord(c, r):
        if r not in replicas:
            replicas[r] = []
        replicas[r].extend(c[3:-3])

    bar = ProgressBar(line_count(filename), msg='Processing trajectory')
    with open(filename) as f:
        current_header = None
        current_coord = []
        for index, line in enumerate(f):
            bar.update(index)
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
    bar.done()

    headers.sort(key=lambda x: x.model)
    headers.sort(key=lambda x: x.replica)
    coordinates = np.array([x for y in sorted(replicas) for x in replicas[y]], np.int16)

    return headers, coordinates

if __name__ == '__main__':
    t, h = read_trajectory('CABS/TRAF', 'CABS/SEQ')
    tra = Trajectory(*t)
    replicas = tra.to_atoms()
    print replicas[0].make_pdb()
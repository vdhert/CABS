"""Module for handling CABS file formats"""

from vector3d import *
import numpy as np
from atom import *
from os.path import join

__all__ = ['CabsLattice']


class CabsLattice:
    """
    This class represent a CABS-like lattice. It is initialized with:
    grid_spacing: distance between grid nodes, default 0.61 A
    r12: tuple with min and max allowed values for CA-CA pseudo-bond length
    r13: tuple with min and max allowed values for CA-CA-CA end distance
    """
    def __init__(self, grid_spacing=0.61, r12=(3.28, 4.27), r13=(4.1, 7.35)):
        self.grid = grid_spacing
        r12min = round((r12[0] / self.grid)**2)
        r12max = round((r12[1] / self.grid)**2)
        r13min = round((r13[0] / self.grid)**2)
        r13max = round((r13[1] / self.grid)**2)
        dim = int(r12max**0.5)

        self.vectors = []
        for i in range(-dim, dim + 1):
            for j in range(-dim, dim + 1):
                for k in range(-dim, dim + 1):
                    l = i * i + j * j + k * k
                    if r12min <= l <= r12max:
                        self.vectors.append(Vector3d(i, j, k))

        n = len(self.vectors)
        self.good = np.zeros((n, n))
        for i in range(n):
            vi = self.vectors[i]
            for j in range(n):
                vj = self.vectors[j]
                if r13min < (vi + vj).mod2() < r13max and vi.cross(vj).mod2():
                    self.good[i, j] = 1

    def cast(self, ch):
        """
        Function that casts a single protein chain onto the lattice.
        Returns a list of tuples with (x, y, z) coordinates of CA atoms.
        """

        if len(ch.atoms) < 3:
            raise ValueError('Protein chain too short!')

        prev = None
        coord = [Vector3d(
            round(ch.atoms[0].coord.x / self.grid),
            round(ch.atoms[0].coord.y / self.grid),
            round(ch.atoms[0].coord.z / self.grid)
        )]

        for atom in ch.atoms[1:]:
            #  iterate over atoms
            min_dr = 1e12
            min_i = -1

            for i, v in enumerate(self.vectors):
                #  iterate over possible vectors

                if len(coord) > 2 and self.good[prev, i] == 0:
                    continue

                new = coord[-1] + v
                dr = (self.grid * new - atom.coord).mod2()
                if dr < min_dr:
                    min_dr = dr
                    min_i = i

            if min_i < 0:
                raise ValueError('Unsolvable geometric problem!')
            else:
                coord.append(coord[-1] + self.vectors[min_i])
                prev = min_i

        coord.insert(0, coord[0] + coord[1] - coord[2])
        coord.append(coord[-1] + coord[-2] - coord[-3])

        return coord

    def make_input_files(self, complex):
        """
        Returns fchains and seq as strings.
        """
        fchains = None
        seq = ''
        for model in complex.models():
            chains = model.chains()
            if not fchains:
                fchains = [''] * len(chains)
                for atom in model:
                    seq += '%5i%1s %1s%3s %1s%3i%6.2f\n' % (
                        atom.resnum,
                        atom.icode,
                        atom.alt,
                        atom.resname,
                        atom.chid,
                        int(atom.occ),
                        atom.bfac
                    )

            for i, chain in enumerate(chains):
                vectors = self.cast(chain)
                fchains[i] += str(len(vectors)) + '\n' + '\n'.join(
                    ['%i %i %i' % (int(v.x), int(v.y), int(v.z)) for v in vectors]
                ) + '\n'
        return ''.join(fchains), seq


if __name__ == '__main__':
    pass
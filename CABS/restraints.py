"""Module for handling distance restraints"""

from CABS.utils import pep2pep1
import random


class Restraint:
    """Class represents single distance restraint"""
    def __init__(self, line, is_side_chain=False):
        i1, i2, d, w = line.split()
        self.id1 = pep2pep1(i1)
        self.id2 = pep2pep1(i2)
        self.distance = float(d)
        self.weight = float(w)
        self.sg = is_side_chain

    def __repr__(self):
        s = '%s %s %.4f %.2f' % (self.id1, self.id2, self.distance, self.weight)
        if self.sg:
            s += ' SG'
        return s

    def update_id(self, ids):
        self.id1 = ids[self.id1]
        self.id2 = ids[self.id2]
        return self


class Restraints:
    """Container for Restraint(s)"""
    def __init__(self, restraints, sg=False):
        self.data = []
        if restraints:
            self.data.extend(Restraint(line, sg) for line in restraints)

    @classmethod
    def from_parser(cls, restraints, sg=False):
        as_strings = [' '.join(str(w) for w in r) for r in restraints]
        return cls(as_strings, sg)

    @classmethod
    def from_file(cls, filename, sg=False):
        with open(filename) as f:
            as_strings = [line for line in f]
        return cls(as_strings, sg)

    def __repr__(self):
        return '\n'.join(str(r) for r in self.data)

    def __iadd__(self, other):
        self.data.extend(other.data)
        return self

    def update_id(self, ids):
        for restr in self.data:
            restr.update_id(ids)
        return self

    def reduce_by(self, factor):
        if 0. < factor < 1.:
            _count = int(len(self.data) * factor)
            self.data = random.sample(self.data, _count)


if __name__ == '__main__':
    pass

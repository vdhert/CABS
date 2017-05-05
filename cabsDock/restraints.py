"""Module for handling distance restraints"""

__all__ = ['Restraint', 'Restraints']


class Restraint:
    """Class represents single distance restraint"""
    def __init__(self, line, is_side_chain=False):
        i1, i2, d, w = line.split()
        self.id1 = i1
        self.id2 = i2
        self.distance = float(d)
        self.weight = float(w)
        self.sg = is_side_chain

    def __repr__(self):
        s = '%s %s %.4f %.2f' % (self.id1, self.id2, self.distance, self.weight)
        if self.sg:
            s += ' SG'
        return s


class Restraints(list):
    """Container for Restraint(s)"""
    def __init__(self, rest, sg=False):
        list.__init__(self)
        if type(rest) is list:
            self.extend([Restraint(line, sg) for line in rest])
        elif type(rest) is str:
            with open(rest) as f:
                self.extend(Restraint(line, sg) for line in f)

    def __repr__(self):
        return '\n'.join(str(r) for r in self)


if __name__ == '__main__':
    pass

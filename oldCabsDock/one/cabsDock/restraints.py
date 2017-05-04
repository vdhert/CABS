"""
Module for handling restraints 
"""


class Restraint:
    def __init__(self, line, sidechains=False):
        i, j, d, w = line.split()
        self.id1 = i
        self.id2 = j
        self.distance = float(d)
        self.weight = float(w)
        self.sg = sidechains


class Restraints:
    def __init__(self, **kwargs):
        self.ca_restraints = []
        self.sg_restraints = []
        if 'from_pdb' in kwargs:
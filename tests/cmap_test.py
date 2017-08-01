from unittest import TestCase
from unittest import main
from cabsDock.utils import SCModeler
import pickle
import numpy
import sys
import random
from argparse import ArgumentParser as AP

class SCModelerTest(TestCase):

    @classmethod
    def setUpClass(cls):
        with open(tdata + 'traj.pck') as f:
            cls.traj = pickle.load(f)
        with open(tdata + 'cplx.pck') as f:
            cls.temp = pickle.load(f)

        cls.aas = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']

        class DummyAtom(object):
            def __init__(self, rn=None):
                self.resname = rn if rn else random.choice(cls.aas)

        cls.dm = DummyAtom

    def test_init(self):
        SCModeler(self.temp.atoms)

    def test_rebuild_one_pickled_CB(self):
        scm = SCModeler(self.temp.atoms)
        frame = self.traj.coordinates[1, 1, ...]
        rbld = scm.rebuild_one(frame, False)
        self.assertEqual(frame.shape, rbld.shape)
        for ca, cb in zip(frame, rbld):
            self.assertLess(abs(numpy.linalg.norm(ca - cb) - 0.7), 0.1)

    def test_rebuild_one_artif_scatter_CB(self):
        for i in self.aas:
            scm = SCModeler([self.dm(), self.dm(i), self.dm()])
            frame = numpy.array([   [0., 0., 0.],
                                    [3.2, 1.5, 0.],
                                    [6.4, 0., 0.]])
            res = scm.rebuild_one(frame, False)
            dist = numpy.linalg.norm(frame[1] - res[1])
            self.assertEqual(0.761, round(dist, 3))

    def test_rebuild_one_artif_compact_CB(self):
        for i in self.aas:
            scm = SCModeler([self.dm(), self.dm(i), self.dm()])
            frame = numpy.array([   [0., 0., 0.],
                                    [3.2, 1.5, 0.],
                                    [5.6, 0., 0.]])
            res = scm.rebuild_one(frame, False)
            dist = numpy.linalg.norm(frame[1] - res[1])
            self.assertEqual(0.740, round(dist, 3))

if __name__ == '__main__':
    ap = AP()
    ap.add_argument('--data_path', default='./tests/data/')
    opt, args = ap.parse_known_args()
    tdata = opt.data_path
    main(argv=sys.argv[:1] + args)
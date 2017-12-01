from unittest import TestCase
from unittest import main
from unittest import SkipTest
from CABS.utils import SCModeler
from CABS.utils import SIDECNT
from CABS.job import FlexTask
from CABS.cmap import ContactMapFactory
from CABS.cmap import ContactMap

import pickle
import numpy
import sys
import os
import random
import tempfile

from argparse import ArgumentParser as AP

class CMapMakerTest(TestCase):

    @classmethod
    def setUpClass(cls):

        class DummyAtom(object):
            n = 0
            def __init__(self, ch):
                self.chid = ch
                self.n = DummyAtom.n
                DummyAtom.n += 1
            def fmt(self):
                return "%s%i" % (self.chid, self.n)

        class DummyAtoms(object):
            atoms = [DummyAtom('A') for i in range(20)]
            atoms += [DummyAtom('B') for i in range(20)]

        class DAs2(object):
            atoms = [DummyAtom('A') for i in range(10)]
            atoms += [DummyAtom('B') for i in range(10)]

        cls.DA = DummyAtoms
        cls.DA2 = DAs2

    def test_factory_init(self):
        cmf = ContactMapFactory('A', 'B', self.DA)
        self.assertEqual(cmf.inds1.__len__(), 20)
        self.assertEqual(cmf.inds2.__len__(), 20)
        self.assertEqual(cmf.ats1.__len__(), 20)
        self.assertEqual(cmf.ats2.__len__(), 20)
        cmf = ContactMapFactory('A', 'AB', self.DA)
        self.assertEqual(cmf.inds1.__len__(), 20)
        self.assertEqual(cmf.inds2.__len__(), 40)
        self.assertEqual(cmf.ats1.__len__(), 20)
        self.assertEqual(cmf.ats2.__len__(), 40)

    def test_factory_mk_dmtx(self):
        cmf1 = ContactMapFactory('A', 'B', self.DA)
        cmf2 = ContactMapFactory('A', 'AB', self.DA)
        cmf3 = ContactMapFactory('A', 'A', self.DA)
        vector = numpy.ones((40, 3)) * numpy.arange(40)[:,numpy.newaxis]
        res1 = cmf1.mk_dmtx(vector)
        res2 = cmf2.mk_dmtx(vector)
        res3 = cmf3.mk_dmtx(vector)
        self.assertEqual(res1.shape, (20, 20))
        self.assertEqual(res2.shape, (20, 40))
        #shapes are okay
        self.assertFalse(len(numpy.diag(res3).nonzero()[0]))
        #zeros on diag in case 3 -- one chain comp. with itself
        self.assertEqual(numpy.max(numpy.diag(res1)), numpy.min(numpy.diag(res1)))
        self.assertEqual(numpy.max(numpy.diag(res1)), numpy.sqrt(20 ** 2 * 3))
        #distance between subsequent res of different artificial "chains"
        #should be equal for all cases and = Sqrt(20 ** 2 * 3)
        bd = numpy.sqrt(3)  # basic distance between artif. res.
        for n, r in enumerate(res3):
            for m, i in enumerate(r):
                self.assertAlmostEqual(i / bd, abs(m - n))
        #subsequent values should be equal to certain number of Sqrt(3)
        #which is dist between subseq. artif. res.

    def test_factory_mk_cmtx(self):
        cmf1 = ContactMapFactory('A', 'B', self.DA)
        for i in range(3, 100):
            arg = numpy.arange(i ** 2)
            mtx = arg.reshape(i, i)
            thr = random.choice(arg)
            res = cmf1.mk_cmtx(mtx, thr)
            self.assertEqual(len(res.nonzero()[0]), thr)
        #TODO: assert Truths on right possitions

    def test_factory_mk_cmap(self):
        frame = numpy.ones((20, 3))
        for i in range(10):
            frame[i,] = frame[i,] + numpy.array([0., i, 0.])
        for i in range(10):
            frame[10 + i,] = frame[10 - i - 1,]
        frame[10:,] = frame[10:,] + numpy.array([0., 0., 1.])
        atraj = numpy.zeros((3, 10, 20, 3))
        for rep in range(3):
            for frm in range(10):
                atraj[rep, frm, :, :] = frame
        for frm, co in enumerate(numpy.linspace(1., 2., 10)):
            atraj[2, frm, :, :] = atraj[2, frm, :, :] * co
        cmf = ContactMapFactory('A', 'B', self.DA2)
        for thr, ncons in ((.9, 0), (1., 0), (1.1, 10), (1.5, 8 * 3 + 2 * 2), (11., 10 ** 2)):
            res = cmf.mk_cmap(atraj, thr)
            self.assertEqual(len(res), 3)
            self.assertTrue(max(map(type, res)) is ContactMap)
            for n, i in enumerate(res[:2]):
                self.assertEqual(i.n, 10)
                self.assertEqual(len(i.cmtx.nonzero()[0]), ncons, "Wrong number of contacts for thr %.2f: %i instead of %i. (Replica %i)" % (thr, len(i.cmtx.nonzero()[0]), ncons, n))
            res1 = cmf.mk_cmap(atraj, thr, replicas=(2,))
            self.assertEqual(res1.__len__(), 1)
            self.assertEqual(res1[0].n, 10)
            ens = len([i for i in numpy.linspace(1., 2., 10) if i < thr])
            # number of frames in which 1A contacts are under the threshold
            if ens > 1:
                self.assertGreaterEqual(len(res1[0].cmtx.nonzero()[0]), ens)
                res2 = cmf.mk_cmap(atraj, thr, replicas=(2,), frames=range(ens))
                self.assertTrue(numpy.array_equal(res1[0].cmtx, res2[0].cmtx))
                res3 = cmf.mk_cmap(atraj, thr, replicas=(2,), frames=range(1, ens))
                self.assertTrue(numpy.array_equal(numpy.clip(res2[0].cmtx - 1, 0, 10), res3[0].cmtx))


class SCModelerTest(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aas = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']

        class DummyAtom(object):
            def __init__(self, rn=None):
                self.resname = rn if rn else random.choice(cls.aas)

        cls.dm = DummyAtom

    def test_init(self):
        SCModeler([self.dm(i) for i in self.aas])
        SCModeler([self.dm(random.choice(self.aas)) for i in range(100)])

    def test_rebuild_one_pickled_CB(self):
        raise Exception('No test')
        #~ scm = SCModeler(self.temp.atoms)
        #~ frame = self.traj.coordinates[1, 1, ...]
        #~ rbld = scm.rebuild_one(frame, False)
        #~ self.assertEqual(frame.shape, rbld.shape)
        #~ for ca, cb in zip(frame, rbld):
            #~ self.assertLess(abs(numpy.linalg.norm(ca - cb) - 0.7), 0.1)

    def mk_random_versor(self, fixed=(None, None, None)):
        v1 = numpy.random.random(3)
        for i, n in zip(fixed, range(3)):
            v1[n] = i if i is not None else v1[n]
        return v1 / numpy.linalg.norm(v1)

    def mk_rot_z_mtx(self, ang=None):
        if not ang:
            rnd = numpy.random.random_integers(2, 8)
            rnd2 = numpy.random.random_integers(1, 2 * rnd -1)
            ang = rnd2 * numpy.pi / rnd
        rot = numpy.array([ [numpy.cos(ang), -1 * numpy.sin(ang), 0],
                            [numpy.sin(ang), numpy.cos(ang), 0],
                            [0, 0, 1]])     #sample rotation
        return rot, ang

    def mk_random_cas(self, scatter=True):
        """Return array of C_alpha vectors in compact or scatter formation."""
        if scatter:
            c2 = [2.65, 2.72, 0.]
            c3 = [5.3, 0., 0]
        else:
            c2 = [3.2, 2.04, 0.]
            c3 = [6.4, 0., 0]
        frame = numpy.array([   [0., 0., 0.],
                                c2,
                                c3])
        vec = self.mk_random_versor() * numpy.random.random_integers(2, 10)
        rot, dummy = self.mk_rot_z_mtx()
        rand_ord = range(3)
        numpy.random.shuffle(rand_ord)
        rot = rot[rand_ord][:,rand_ord]
        return numpy.array([numpy.dot(rot, i) + vec for i in frame]), vec, rot

    def test_mk_local_system(self):
        mtx = numpy.array([ [0., 0., 0.],
                            [3.2, 1.5, 0.],
                            [6.4, 0., 0.]]) # sample CA positions (scatter)
        ang = numpy.pi / 6
        rot, ang = self.mk_rot_z_mtx(ang)
        c1, c2, c3 = [numpy.dot(rot, i) for i in mtx]
        x, y, z, dst = SCModeler._mk_local_system(c1, c2, c3)
        self.assertEqual(dst, 6.4)
        for i in (x, y, z):
            self.assertEqual(1, numpy.linalg.norm(i))   # we expect versors
            for j in (x, y, z):
                if i is j: continue
                self.assertAlmostEqual(numpy.dot(i, j), 0)  # versors should be orthogonal

    def test_calc_nodes_line(self):
        for i in range(100):
            v1 = self.mk_random_versor()
            rot, ang = self.mk_rot_z_mtx()
            v2 = numpy.dot(rot, v1)
            nv = SCModeler._calc_nodes_line(v1, v2)
            self.assertAlmostEqual(1, numpy.linalg.norm(nv))
            self.assertAlmostEqual(0, numpy.dot(nv, v1))
            self.assertAlmostEqual(0, numpy.dot(nv, v2))

    def test_calc_trig_fnc(self):
        for i in range(100):
            v1 = self.mk_random_versor((None, None, 0))
            rot, ang = self.mk_rot_z_mtx()
            v2 = numpy.dot(rot, v1)
            cos, sin = SCModeler._calc_trig_fnc(v1, v2, numpy.array((0, 0, 1)))
            xx = numpy.linalg.norm(numpy.cross(v1, v2))
            self.assertAlmostEqual(sin, numpy.sin(ang))
            self.assertAlmostEqual(cos, numpy.cos(ang))

    def test_calc_rot_mtx(self):
        for i in range(100):
            scm = SCModeler([self.dm(), self.dm(), self.dm()])
            frame, vec, rot = self.mk_random_cas()
            rot, dist = scm._calc_rot_mtx(*frame)
            self.assertEqual(numpy.linalg.norm(frame[0] - frame[2]), dist)
            xp = numpy.dot(rot, (1, 0, 0))
            yp = numpy.dot(rot, (0, 1, 0))
            zp = numpy.dot(rot, (0, 0, 1))
            for tv in (xp, yp, zp):
                self.assertAlmostEqual(numpy.linalg.norm(tv), 1)
            for v in SIDECNT.values():
                v1 = v[:3]
                v2 = v[3:]
                n1 = numpy.linalg.norm(v1)
                n2 = numpy.linalg.norm(v2)
                n1p = numpy.linalg.norm(numpy.dot(rot, v1))
                n2p = numpy.linalg.norm(numpy.dot(rot, v2))
                self.assertAlmostEqual(n1, n1p)
                self.assertAlmostEqual(n2, n2p)

    def test_coef(self):
        for k, v in SIDECNT.items():
            if k == 'GLY': continue
            for i in numpy.arange(5, 7, .1):
                res = SCModeler._calc_scatter_coef(i)
                self.assertGreaterEqual(res, 0)
                self.assertLessEqual(res, 1)
                v1 = numpy.array(v[:3])
                v2 = numpy.array(v[3:])
                diff = numpy.linalg.norm(v1 - v2)
                av_vec = (v1 * res + v2 * (1 - res))
                av_diff1 = numpy.linalg.norm(av_vec - v1)
                av_diff2 = numpy.linalg.norm(av_vec - v2)
                self.assertLessEqual(av_diff1, diff)
                self.assertLessEqual(av_diff2, diff)

    def test_rebuild_one_artif_scatter_CB(self):
        for dmm in range(50):
            for i in self.aas:
                scm = SCModeler([self.dm(), self.dm(i), self.dm()])
                frame, dummy1, dummy2 = self.mk_random_cas()
                res = scm.rebuild_one(frame, False)
                dist = numpy.linalg.norm(frame[1] - res[1])
                prop_dst = numpy.linalg.norm(SIDECNT['ALA'][3:])
                #~ if round(prop_dst, 3) != round(dist, 3):
                    #~ scm._calc_rot_mtx(*frame, dbg=True)

                self.assertAlmostEqual(round(prop_dst, 3), round(dist, 3))

    def test_rebuild_one_artif_compact_CB(self):
        for dmm in range(50):
            for i in self.aas:
                scm = SCModeler([self.dm(), self.dm(i), self.dm()])
                frame, dmm1, dmm2 = self.mk_random_cas()
                res = scm.rebuild_one(frame, False)
                dist = numpy.linalg.norm(frame[1] - res[1])
                prop_dst = numpy.linalg.norm(SIDECNT['ALA'][:3])
                self.assertEqual(round(prop_dst, 3), round(dist, 3))

if __name__ == '__main__':
    ap = AP()
    ap.add_argument('--data_path', default='./tests/data/')
    ap.add_argument('--fast', default=False, action='store_true')
    cla, args = ap.parse_known_args()
    main(argv=sys.argv[:1] + args)

from unittest import TestCase
from unittest import main
from unittest import SkipTest
from cabsDock.utils import SCModeler
from cabsDock.utils import SIDECNT
from cabsDock.job import Job
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
        if cla.fast:
            raise SkipTest
        ddir = tempfile.gettempdir() + "/tmpCABS/"
        try: os.mkdir(ddir)
        except OSError: pass
        cls.j = Job(
            work_dir=ddir,
            receptor='1klu:AB',
            ligand=[['GELIGTLNAAKVPAD:CCCEEEECCEECCCC', 'random', 'random']],
            mc_cycles=20,
            mc_steps=50,
            replicas=10,
            load_cabs_files=(cla.data_path + '1klu/TRAF', cla.data_path + '1klu/SEQ'),
            AA_rebuild=False
            )
        cls.j.cabsdock()

    def test_contact_frequencies(self):
        pass

class SCModelerTest(TestCase):

    @classmethod
    def setUpClass(cls):
        with open(cla.data_path + 'traj.pck') as f:
            cls.traj = pickle.load(f)
        with open(cla.data_path + 'cplx.pck') as f:
            cls.temp = pickle.load(f)

        cls.aas = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']

        class DummyAtom(object):
            def __init__(self, rn=None):
                self.resname = rn if rn else random.choice(cls.aas)

        cls.dm = DummyAtom

    def test_init(self):
        SCModeler(self.temp.atoms)

    #~ def test_rebuild_one_pickled_CB(self):
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
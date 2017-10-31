from unittest import TestCase
from unittest import main
from unittest import SkipTest
import tempfile
import os
import sys
from argparse import ArgumentParser as AP
from shutil import rmtree

from CABS.job import FlexTask
from CABS.optparser import FlexParser

ap = AP()
ap.add_argument('--fast', action='store_true', default=False)


class FlexTest(TestCase):

    def mk_parser(self, lst):
        ddir = tempfile.gettempdir() + "/flexInitCABS"
        prsr = FlexParser.parse_args(lst + ['--work-dir', ddir])
        return ddir, prsr

    def test_init(self):
        ddir, prsr = self.mk_parser(['-i', '1hpw', '--verbose', '-1'])
        tsk = FlexTask(**vars(prsr))
        rmtree(ddir)

    def test_default(self):
        ddir, prsr = self.mk_parser(['-i', '1hpw', '--verbose', '-1'])
        tsk = FlexTask(**vars(prsr))
        self.assertEqual(tsk.mc_steps, 50)
        self.assertEqual(tsk.replicas, 1)
        self.assertEqual(tsk.mc_annealing, 20)
        self.assertEqual(tsk.mc_cycles, 50)
        self.assertEqual(tsk.input_protein, '1hpw')
        self.assertEqual(tsk.protein_restraints, ('ss2', 3, 3.8, 8.0))
        self.assertEqual(tsk.work_dir, ddir)
        self.assertEqual(tsk.verbose, -1)
        rmtree(ddir)
        #TODO test no input given

    def test_output_completeness(self):
        ddir, prsr = self.mk_parser([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--verbose', '-1',
        ])
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        lst = os.listdir(ddir)
        outSet = {  'output_pdbs':
                        ['%s_%i.pdb' % (tp, i) for tp in ('model', 'cluster') for i in range(10)] + ['replica.pdb'],
                    'output_data':
                        ['all_rmsds_A.txt', 'lowest_rmsds_A.txt', 'medoids_rmsds_A.txt', 'reference_alignment_target.csv'],
                    'plots':
                        ['%s.%s' % (n, e) for e in ('csv', 'svg') for n in ('E_RMSD_A_total', 'RMSD_frame_A_replica_0')]
                    }
        self.assertEqual(set(lst).intersection(outSet), set(outSet.keys()))
        for i, j in outSet.items():
            lst = os.listdir(ddir + '/' + i)
            self.assertEqual(set(lst).intersection(j), set(j))
        rmtree(ddir)

    def test_run_default(self):
        if cla.fast:
            raise SkipTest
        ddir, prsr = self.mk_parser(['-i', '1hpw', '--verbose', '-1'])
        tsk = FlexTask(**vars(prsr))
        tsk.run()

        def clc_len(fn):
            with open(fn) as f:
                return len(f.raedlines())

        csvs = {i: clc_len(ddir + '/plots/' + i) for i in os.listdir(ddir + '/plots') if 'csv' in i}
        for i, j in csvs.items():
            self.assertEqual(j, 129 if i == 'RMSF.csv' else 1000)

        rmtree(ddir)

    def test_run(self):
        ddir, prsr = self.mk_parser([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--verbose', '-1',
        ])
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        rmtree(ddir)

    def test_save_CABS_files(self):
        ddir, prsr = self.mk_parser([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--save-cabs-files',
            '--verbose', '-1',
        ])
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        self.assertTrue('SEQ' in os.listdir(ddir))
        self.assertTrue('TRAF' in os.listdir(ddir))
        rmtree(ddir)

    def test_analysis(self):
        ddir, prsr = self.mk_parser([
            '-i', '1hpw',
            '--load-cabs-files', 'tests/data/1hpw/TRAF', 'tests/data/1hpw/SEQ',
            '--verbose', '-1',
        ])
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        rmtree(ddir)
        #TODO add files to data/1hpw
        #TODO assert results are the same as in data


if __name__ == '__main__':
    cla, args = ap.parse_known_args()
    main(argv=sys.argv[:1] + args)
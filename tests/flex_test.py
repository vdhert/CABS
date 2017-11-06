from unittest import TestCase
from unittest import main
from unittest import SkipTest
import tempfile
import os
import sys
from argparse import ArgumentParser as AP
from shutil import rmtree
from StringIO import StringIO

from CABS.job import FlexTask
from CABS.optparser import FlexParser

ap = AP()
ap.add_argument('--fast', action='store_true', default=False)


def add_args_and_workdir(args):
    def decorTD(mth):
        def wrapped(self):
            ddir, prsr = self.mk_parser(args)
            mth(self, ddir, prsr)
            rmtree(ddir)
        return wrapped
    return decorTD

class FlexParserTest(TestCase):

    def test_FPinit(self):
        #~ FlexParser.parse_args([])
        FlexParser.parse_args(['-i', '1hpw'])
        #~ std = StringIO()
        #~ sys.stdout = std
        #~ FlexParser.parse_args(['-i',])
        #~ sys.stdout = sys.__stdout__
        #~ self.assertTrue('error' in std.read())
        #TODO -- przechwytywanie stdout i sprawdzanie wyjatku przy braku inputu


class FlexTest(TestCase):

    def mk_parser(self, lst):
        ddir = tempfile.gettempdir() + "/flexInitCABS"
        prsr = FlexParser.parse_args(lst + ['--work-dir', ddir])
        return ddir, prsr

    @add_args_and_workdir(['-i', '1hpw', '--verbose', '-1'])
    def test_FTinit(self, ddir, prsr):
        tsk = FlexTask(**vars(prsr))

    @add_args_and_workdir(['-i', '1hpw', '--verbose', '-1'])
    def test_default(self, ddir, prsr):
        tsk = FlexTask(**vars(prsr))
        self.assertEqual(tsk.mc_steps, 50)
        self.assertEqual(tsk.replicas, 1)
        self.assertEqual(tsk.mc_annealing, 20)
        self.assertEqual(tsk.mc_cycles, 50)
        self.assertEqual(tsk.input_protein, '1hpw')
        self.assertEqual(tsk.protein_restraints, ('ss2', 3, 3.8, 8.0))
        self.assertEqual(tsk.work_dir, ddir)
        self.assertEqual(tsk.verbose, -1)

    @add_args_and_workdir([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--verbose', '-1',
        ])
    def test_output_completeness(self, ddir, prsr):
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

    @add_args_and_workdir(['-i', '1hpw', '--verbose', '-1'])
    def test_run_default(self, ddir, prsr):
        if cla.fast:
            raise SkipTest
        tsk = FlexTask(**vars(prsr))
        tsk.run()

        def clc_len(fn):
            with open(fn) as f:
                return len(f.raedlines())

        csvs = {i: clc_len(ddir + '/plots/' + i) for i in os.listdir(ddir + '/plots') if 'csv' in i}
        for i, j in csvs.items():
            self.assertEqual(j, 129 if i == 'RMSF.csv' else 1000)

        rmtree(ddir)

    @add_args_and_workdir([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--verbose', '-1',
        ])
    def test_run(self, ddir, prsr):
        tsk = FlexTask(**vars(prsr))
        tsk.run()

    @add_args_and_workdir([
            '-i', '1hpw',
            '--mc-steps', '1',
            '--mc-cycles', '10',
            '--mc-annealing', '1',
            '--save-cabs-files',
            '--verbose', '-1',
        ])
    def test_save_CABS_files(self, ddir, prsr):
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        self.assertTrue('SEQ' in os.listdir(ddir))
        self.assertTrue('TRAF' in os.listdir(ddir))

    @add_args_and_workdir([
            '-i', '1hpw',
            '--load-cabs-files', 'tests/data/1hpw/TRAF', 'tests/data/1hpw/SEQ',
            '--verbose', '-1',
        ])
    def test_analysis(self, ddir, prsr):
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        with open(os.path.join(ddir, 'output_pdbs', 'replica.pdb')) as f:
            result = f.readlines()
        with open('tests/data/1hpw/replica.pdb') as f:
            template = f.readlines()
        self.assertEqual(result, template)
        #TODO jak bedzie seed -- dorobic klastrowanie
        #TODO zbadac inne pliki -- jakie?

    @add_args_and_workdir([
            '-i', '2gb1',
            '--verbose', '-1',
        ])
    def test_rmsf(self, ddir, prsr):
        """Plik tests/data/2gb1.rmsf zawiera znormalizowane do max rmsf z wersji serwerowej."""
        tsk = FlexTask(**vars(prsr))
        tsk.run()
        with open(os.path.join(ddir, 'plots/RMSF.csv')) as f:
            rmsf = [float(l.split()[1]) for l in f.readlines()]
        rmsf = [i/max(rmsf) for i in rmsf]
        with open('tests/data/2gb1.rmsf') as f:
            ref = [float(l.split()[1]) for l in f.readlines()]
        diffs = [i - j for i, j in zip(rmsf, ref)]
        av = sum(map(abs, diffs)) / len(diffs)
        self.assertLess(av, .15)    # average lower than .15
        mavs = [sum(diffs[i: i+6])/6. for i in range(len(diffs) - 6)]
        self.assertGreater(len([i for i in map(abs, mavs) if i < .2]), .8 * len(diffs))
        # more than 80 % of residues 6-moving-avs is below .2


if __name__ == '__main__':
    cla, args = ap.parse_known_args()
    main(argv=sys.argv[:1] + args)
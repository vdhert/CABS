from cabsDock.job import Job
from matplotlib import pyplot
import numpy as np
from os.path import join
import copy_reg
import types


def _reduce_method(meth):
    return getattr, (meth.__self__, meth.__func__.__name__)


copy_reg.pickle(types.MethodType, _reduce_method)


class Benchmark(object):
    """
    Class for testing the performance against a reference set of structures (-> cases)
    """

    def __init__(self, benchmark_file=None):
        """

        :param benchmark_file: text file with lines in format: bound_pdb_code, bound_receptor_chain_id,
        bound_peptide_chain_id, peptide_sequence, peptide_secondary_structure,
        """
        super(Benchmark, self).__init__()
        self.benchmark_rmsds_10 = []
        self.benchmark_rmsds_1k = []
        self.benchmark_rmsds_10k = []
        self.cases = []

        self.benchmark_file = benchmark_file
        f = open(self.benchmark_file, 'r')
        lines = f.readlines()

        for line in lines:
            row = line.split()
            self.cases.append(
                Case(
                    bound_pdb_code=row[0],
                    bound_receptor_chain_id=row[1],
                    bound_peptide_chain_id=row[2],
                    peptide_sequence=row[3],
                    peptide_secondary_structure=row[4],
                    unbound_pdb_code=None,
                    unbound_receptor_chain_id=None,
                    unbound_peptide_chain_id=None
                )
            )

    def __iter__(self):
        return self.cases.__iter__()

    def bench_set(self):
        for element in self.cases:
            element.setup_case(test=True)

    def bench_run(self, parallel=False):
        if parallel:
            pass
        else:
            for element in self.cases:
                if element.is_valid:
                    element.run_case()

    def bench_analyze(self, histograms=False):
        self.benchmark_lowest_rmsds_10k = []
        self.benchmark_lowest_rmsds_1k = []
        self.benchmark_lowest_rmsds_10 = []
        for element in self.cases:
            rmsds = element.rmsds
            if histograms:
                pyplot.hist(rmsds[0])
                pyplot.hist(rmsds[1])
                pyplot.hist(rmsds[2])
                pyplot.savefig(str(element.bound_pdb_code) + '.png')
                pyplot.clf()
            self.benchmark_rmsds_10k += rmsds[0]
            self.benchmark_rmsds_1k += rmsds[1]
            self.benchmark_rmsds_10 += rmsds[2]
            self.benchmark_lowest_rmsds_10k.append(rmsds[3])
            self.benchmark_lowest_rmsds_1k.append(rmsds[4])
            self.benchmark_lowest_rmsds_10.append(rmsds[5])
        # print the performance statistics on this benchmark
        number_of_cases = float(len([True for case in self.cases if case.is_valid]))
        high_quality_10k = len([rms for rms in self.benchmark_lowest_rmsds_10k if rms < 3.0])
        med_quality_10k = len([rms for rms in self.benchmark_lowest_rmsds_10k if 3.0 <= rms <=
                               5.5])
        low_quality_10k = len([rms for rms in self.benchmark_lowest_rmsds_10k if rms > 5.5])

        print("1k")
        high_quality_1k = len([rms for rms in self.benchmark_lowest_rmsds_1k if rms < 3.0])
        med_quality_1k = len([rms for rms in self.benchmark_lowest_rmsds_1k if 3.0 <= rms <=
                              5.5])
        low_quality_1k = len([rms for rms in self.benchmark_lowest_rmsds_1k if rms > 5.5])
        print("top10")
        high_quality_10 = len([rms for rms in self.benchmark_lowest_rmsds_10 if rms < 3.0])
        med_quality_10 = len([rms for rms in self.benchmark_lowest_rmsds_10 if 3.0 <= rms <=
                              5.5])
        low_quality_10 = len([rms for rms in self.benchmark_lowest_rmsds_10 if rms > 5.5])
        if histograms:
            pyplot.hist(self.benchmark_rmsds_10k)
            pyplot.hist(self.benchmark_rmsds_1k)
            pyplot.hist(self.benchmark_rmsds_10)
            pyplot.savefig('benchmark.png')
            pyplot.clf()
        print(
            'Statistics for 10k\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_10k,
                high_quality_10k / number_of_cases,
                med_quality_10k,
                med_quality_10k / number_of_cases,
                low_quality_10k,
                low_quality_10k / number_of_cases
            )
        )
        print(
            'Statistics for 1k\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_1k,
                high_quality_1k / number_of_cases,
                med_quality_1k,
                med_quality_1k / number_of_cases,
                low_quality_1k,
                low_quality_1k / number_of_cases
            )
        )
        print(
            'Statistics for 10\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_10,
                high_quality_10 / number_of_cases,
                med_quality_10,
                med_quality_10 / number_of_cases,
                low_quality_10,
                low_quality_10 / number_of_cases
            )
        )

    def multi_run(self, histograms=False):
        import multiprocessing
        import sys
        import os

        sys.path.append(os.getcwd())
        print('Initiating pool. Cpu count is {0}'.format(multiprocessing.cpu_count()))
        pool = multiprocessing.Pool()
        output = {}
        for element in self.cases:
            if element.is_valid:
                output[element.bound_pdb_code] = pool.apply_async(element.run_case)
        pool.close()
        pool.join()
        self.benchmark_rmsds_10k = []
        self.benchmark_rmsds_1k = []
        self.benchmark_rmsds_10 = []
        self.benchmark_lowest_rmsds_10k = []
        self.benchmark_lowest_rmsds_1k = []
        self.benchmark_lowest_rmsds_10 = []
        for bialko in output.keys():
            rmsds, work_dir = output[bialko].get()
            if histograms:
                pyplot.hist(rmsds[0])
                pyplot.hist(rmsds[1])
                pyplot.hist(rmsds[2])
                pyplot.savefig(str(bialko) + '.png')
                pyplot.clf()
            self.benchmark_rmsds_10k += rmsds[0]
            self.benchmark_rmsds_1k += rmsds[1]
            self.benchmark_rmsds_10 += rmsds[2]
            self.benchmark_lowest_rmsds_10k.append(rmsds[3])
            self.benchmark_lowest_rmsds_1k.append(rmsds[4])
            self.benchmark_lowest_rmsds_10.append(rmsds[5])
            np.save(join(work_dir, 'rmsds'), rmsds)
            f = open(join(work_dir, 'best_rmsd.txt'), 'w')
            f.write(
                '{0} {1} {2} \n'.format(rmsds[3], rmsds[4], rmsds[5])
            )
            f.close()

        number_of_cases = float(
            len(self.cases)
        )
        high_quality_10k = len(
            [rms for rms in self.benchmark_lowest_rmsds_10k if rms < 3.0]
        )
        med_quality_10k = len(
            [rms for rms in self.benchmark_lowest_rmsds_10k if 3.0 <= rms <= 5.5]
        )
        low_quality_10k = len(
            [rms for rms in self.benchmark_lowest_rmsds_10k if rms > 5.5]
        )

        high_quality_1k = len(
            [rms for rms in self.benchmark_lowest_rmsds_1k if rms < 3.0]
        )
        med_quality_1k = len(
            [rms for rms in self.benchmark_lowest_rmsds_1k if 3.0 <= rms <= 5.5]
        )
        low_quality_1k = len(
            [rms for rms in self.benchmark_lowest_rmsds_1k if rms > 5.5]
        )

        high_quality_10 = len(
            [rms for rms in self.benchmark_lowest_rmsds_10 if rms < 3.0]
        )
        med_quality_10 = len(
            [rms for rms in self.benchmark_lowest_rmsds_10 if 3.0 <= rms <= 5.5]
        )
        low_quality_10 = len(
            [rms for rms in self.benchmark_lowest_rmsds_10 if rms > 5.5]
        )
        if histograms:
            pyplot.hist(self.benchmark_rmsds_10k)
            pyplot.hist(self.benchmark_rmsds_1k)
            pyplot.hist(self.benchmark_rmsds_10)
            pyplot.savefig('benchmark.png')
            pyplot.clf()
        print(
            'Statistics for 10k\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_10k,
                100 * (high_quality_10k / number_of_cases),
                med_quality_10k,
                100 * (med_quality_10k / number_of_cases),
                low_quality_10k,
                100 * (low_quality_10k / number_of_cases)
            )
        )
        print(
            'Statistics for 1k\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_1k,
                100 * (high_quality_1k / number_of_cases),
                med_quality_1k,
                100 * (med_quality_1k / number_of_cases),
                low_quality_1k,
                100 * (low_quality_1k / number_of_cases)
            )
        )
        print(
            'Statistics for 10\nhigh quality: {0} ({1}%) \nmedium quality: {2} ({3}%) \nlow quality: {4} ({5}%)'.format(
                high_quality_10,
                100 * (high_quality_10 / number_of_cases),
                med_quality_10,
                100 * (med_quality_10 / number_of_cases),
                low_quality_10,
                100 * (low_quality_10 / number_of_cases)
            )
        )


class Case(object):
    """
    Case represents a run and the reference structure it should be compared to.
    """

    def __init__(
            self,
            bound_pdb_code=None,
            bound_receptor_chain_id=None,
            bound_peptide_chain_id=None,
            peptide_sequence=None,
            peptide_secondary_structure=None,
            unbound_pdb_code=None,
            unbound_receptor_chain_id=None,
            unbound_peptide_chain_id=None
    ):
        super(Case, self).__init__()
        self.bound_pdb_code = bound_pdb_code
        self.bound_receptor_chain_id = bound_receptor_chain_id
        self.bound_peptide_chain_id = bound_peptide_chain_id
        self.peptide_sequence = peptide_sequence
        self.peptide_secondary_structure = peptide_secondary_structure
        self.unbound_pdb_code = unbound_pdb_code
        self.unbound_receptor_chain_id = unbound_receptor_chain_id
        self.unbound_peptide_chain_id = unbound_peptide_chain_id
        self.rmsds = []
        self.is_valid = False

    def setup_case(self, test):
        print('Processing job {0}'.format(self.bound_pdb_code))
        try:
            if test:
                self.job = Job(receptor=self.bound_pdb_code + ':' + self.bound_receptor_chain_id,
                               ligand=[[self.peptide_sequence]],
                               work_dir=self.bound_pdb_code,
                               mc_cycles=1,
                               replicas=10,
                               mc_steps=1,
                               native_pdb=self.bound_pdb_code,
                               native_receptor_chain=self.bound_receptor_chain_id,
                               native_peptide_chain=self.bound_peptide_chain_id,
                               )
            else:
                self.job = Job(receptor=self.bound_pdb_code + ':' + self.bound_receptor_chain_id,
                               ligand=[[self.peptide_sequence]],
                               work_dir=self.bound_pdb_code,
                               mc_cycles=50,
                               mc_steps=50,
                               replicas=10,
                               native_pdb=self.bound_pdb_code,
                               native_receptor_chain=self.bound_receptor_chain_id,
                               native_peptide_chain=self.bound_peptide_chain_id,
                               )
        except Exception as errr:
            print(
                '... NOT DONE.\nJob {0} resulted in an ERROR. The ERROR message:\n{1}'
                    .format(self.bound_pdb_code, errr.message)
            )
        else:
            self.is_valid = True
            print('... done.')

    def run_case(self):
        if not self.is_valid:
            raise Exception('The {0} job is not valid for run.'.format(self.bound_pdb_code))
        results = self.job.run_job()
        self.rmsds = [
            results['rmsds_10k'],
            results['rmsds_1k'],
            results['rmsds_10'],
            results['lowest_10k'],
            results['lowest_1k'],
            results['lowest_10']
        ]
        self.work_dir = self.job.config['work_dir']
        print(self.bound_pdb_code, self.rmsds[5])
        return self.rmsds, self.work_dir


from cabsDock.job import Job
from utils import next_letter
import numpy as np
from os.path import join


class Benchmark(object):
    """ docstring for Benchmark """

    def __init__(self, benchmark_file=None):
        super(Benchmark, self).__init__()
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
                element.run_case()


class Case(object):
    """ docstring for Case """

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

    def setup_case(self, test=True):
        print('Processing job {0}'.format(self.bound_pdb_code))
        try:
            if test:
                self.job = Job(receptor=self.bound_pdb_code + ':' + self.bound_receptor_chain_id,
                               ligand=[[self.peptide_sequence]],
                               work_dir=self.bound_pdb_code,
                               mc_cycles=1,
                               replicas=1,
                               mc_steps=1,
                               native_pdb=self.bound_pdb_code,
                               native_receptor_chain=self.bound_receptor_chain_id,
                               native_peptide_chain=self.bound_peptide_chain_id,
                               model_peptide_chain=next_letter(self.bound_receptor_chain_id)
                               )
            else:
                # noinspection PyAttributeOutsideInit
                self.job = Job(receptor=self.bound_pdb_code + ':' + self.bound_receptor_chain_id,
                               ligand=[[self.peptide_sequence]],
                               work_dir=self.bound_pdb_code,
                               mc_cycles=1,
                               mc_steps=5,
                               replicas=10,
                               native_pdb=self.bound_pdb_code,
                               native_receptor_chain=self.bound_receptor_chain_id,
                               native_peptide_chain=self.bound_peptide_chain_id,
                               model_peptide_chain=next_letter(self.bound_receptor_chain_id)
                               )
        except Exception as errr:
            print(
                '... NOT DONE.\nJob {0} resulted in an ERROR. The ERROR message:\n{1}'
                .format(self.bound_pdb_code, errr.message)
            )
        else:
            print('... done.')

    def run_case(self):
        rmsds_10k, rmsds_1k, rmsds_10, work_dir = self.job.run_job()
        np.save(join(work_dir, 'rmsds_10k'), rmsds_10)
        np.save(join(work_dir, 'rmsds_1k'), rmsds_1k)
        np.save(join(work_dir, 'rmsds_10'), rmsds_10)
        f = open(join(work_dir, 'best_rmsd.txt'), 'w')
        f.write(
            str(sorted(rmsds_10)[0]) + ' ' + str(sorted(rmsds_1k)[0]) + ' ' + str(sorted(rmsds_10k)[0]) + '\n'
        )
        f.close()


bench = Benchmark(benchmark_file="/Users/maciek/Desktop/lista_kompl.txt")
bench.cases = bench.cases[0:3]
print([case.bound_pdb_code for case in bench.cases])
bench.bench_set()
bench.bench_run()

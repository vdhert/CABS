#!/bin/env python

#SBATCH --job-name=multiprocess
#SBATCH --nodes=3
#SBATCH --exclusive
#SBATCH --ntasks-per-node=10

import sys
import os

from cabsDock.benchmark import Benchmark
sys.path.append(os.getcwd())
bench = Benchmark(benchmark_file="./lista_kompl.txt")
bench.bench_set()
bench.cases = bench.cases[14:16]
print([case.bound_pdb_code for case in bench.cases])
bench.bench_run()
bench.bench_analyze()
#bench.multi_run()

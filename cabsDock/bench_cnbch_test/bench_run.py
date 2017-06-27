#!/bin/env python

#SBATCH --job-name=bench
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=10

import sys
import os

from cabsDock.benchmark import Benchmark
sys.path.append(os.getcwd())
bench = Benchmark(benchmark_file="./lista_cut.txt")
bench.bench_set(test=True)
print([case.bound_pdb_code for case in bench.cases])
bench.cases = bench.cases[1:2]
bench.bench_run()
bench.bench_analyze()
# bench.multi_run()
# bench.bench_analyze()

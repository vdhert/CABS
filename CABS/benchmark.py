from glob import glob
from os import mkdir, getcwd
from shutil import copyfile


import time
from subprocess import call

from CABS.pbsgen import PbsGenerator


class MBRunner(object):
    def run_standard(self):
        with open('benchmarking_logfile_{}.txt'.format(time.strftime('%x').replace('/', '')), 'a+b') as log:
            print("Bound benchmarks:")
            log.write("MB Bound benchmarks:\n")
            for i in range(1,4):
                br = BenchmarkRunner(benchmark_file='./benchmark_data/MB_bench_bound_{}.txt'.format(i), runtype='bound')
                time.sleep(1)
                command = br.run_benchmark(test=True)
                log.write(command)
                log.write(';\n')
            print("Unbound benchmarks:")
            log.write("MB Unbound benchmarks:\n")
            for i in range(1,4):
                br = BenchmarkRunner(benchmark_file='./benchmark_data/MB_bench_unbound_{}.txt'.format(i), runtype='unbound')
                time.sleep(1)
                command = br.run_benchmark(test=True)
                log.write(command)
                log.write(';\n')

class StandardRunner(object):
    def run_standard(self):
        with open('benchmarking_logfile_{}.txt'.format(time.strftime('%x').replace('/', '')), 'a+b') as log:
            print("Bound benchmarks:")
            log.write("Standard Bound benchmarks:\n")
            for i in xrange(3):
                br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_bound_cases.txt', runtype='bound')
                time.sleep(1)
                command = br.run_benchmark(test=True)
                log.write(command)
                log.write(';\n')
            print("Unbound benchmarks:")
            log.write("Standard Unbound benchmarks:\n")
            for i in xrange(3):
                br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_unbound_cases.txt', runtype='unbound')
                time.sleep(1)
                command = br.run_benchmark(test=True)
                log.write(command)
                log.write(';\n')

    def run_standard_flex(self, version = 'multi', name_prefix=''):
        with open('flex_benchmarking_logfile_{}.txt'.format(time.strftime('%x').replace('/', '')), 'a+b') as log:
            br = BenchmarkRunner(benchmark_file='./benchmark_data/cabsflex_'+version+'chain.txt', mode='cabsflex', runtype='flex', name_prefix=name_prefix)
            command = br.run_benchmark(test=True)
            log.write(command)
            log.write(';\n')

    def __init__(self):
        pass


class BenchmarkRunner(object):
    def __init__(self, benchmark_file, options={'-v': '4', '--remote':''}, name='', name_prefix='', runtype='bound', mode='cabsdock'):
        self.benchmark_file = benchmark_file
        self.options = options
        self.name = name
        self.runtype = runtype
        self.mode = mode
        self.name_prefix = name_prefix

    def setup(self):
        if self.name == '':
            benchdir = '/STORAGE/DATA/mciem'+'/'+self.name_prefix+'benchrun_'
        else:
            benchdir = self.name
        benchdir+=time.strftime("%c").replace(' ','_')+''
        self.benchdir = benchdir
        try:
            mkdir(benchdir)
        except OSError:
            pass
        self.pbsgen = PbsGenerator(benchmark_list=self.benchmark_file, nonstandard_options_dict=self.options,
                                   rundir=self.benchdir+'/run', runtype=self.runtype)
        self.pbsgen.pbs_script(scriptdir=benchdir+'/pbs')


    def save_log(self):
        with open(self.benchdir+'/logfile', 'w') as logfile:
            logfile.write('#Run started: '+time.strftime("%c"))
            options_as_str = '\n'.join([str(key)+' = '+str(val) for (key,val) in self.options.items()])+'\n' if self.options else ''
            logfile.write('\n#Used options:\n'+options_as_str)
            cases_as_str = '\n'.join([case.work_dir for case in self.pbsgen.cases])
            logfile.write('#Cases:\n'+cases_as_str)

    def run_benchmark(self, test=False, qsub_options={'-l walltime':'60:00:00'}):
        self.setup()
        self.save_log()
        options_as_str = ' '.join([str(key)+'='+str(val) for (key, val) in qsub_options.items()])
        command = 'for f in {}/pbs/*.pbs; do qsub {} $f; done'.format(self.benchdir, options_as_str)
        print command
        if test:
           return command
        else:
            call(command)

class BenchmarkAnalyser(object):
    def __init__(self, done_benchdir):
        self.done_benchdir = done_benchdir
        self.cases = []
        self.read_log()


    def read_log(self):
        benchlog = glob(self.done_benchdir+'/logfile')
        if len(benchlog) != 1:
            print benchlog
            raise Exception('Broken or missing (or multiple!) logfile.')
        else:
            self.benchlog = benchlog[0]
        print benchlog
        with open(self.benchlog) as log:
            writecase = False
            for line in log:
                if not line.strip():
                    print "Empty line"
                else:
                    print line
                    thisline = line.split()
                    if 'Cases' in thisline[0]:
                        writecase = True
                    elif writecase and '#' in thisline[0]:
                        writecase = False
                    elif writecase:
                        self.cases.append(thisline[0].upper())

    def get_pearsons(self):
        for case in self.cases:
            print glob('./benchmark_data/cabsflex_reference_rmsf/1A6W_*')
            try:
                reffile = open('./benchmark_data/cabsflex_reference_rmsf/'+case+'_rmsf.txt')
                reference_rmsf = [line for line in reffile.readlines()]
                print reference_rmsf
            except OSError:
                pass

    def read_rmsds(self):
        benchmark_results = {
            'lowest_10k': [],
            'lowest_1k': [],
            'lowest_10': [],
        }
        self.successful_runs = 0.
        with open(self.done_benchdir+'/lowest_rmsds_benchmark.txt', 'w') as outfile:
            for case in self.cases:
                lowest_rmsd_files = glob(self.done_benchdir + '/run/' + case + '/output_data/lowest_rmsds*')
                try:
                    lowest_rmsds = open(lowest_rmsd_files[0])
                except:
                    print('invalid rmsd file for case {}'.format(case))
                    outfile.write(';'.join([case, '---', '---', '---'])+'\n')
                    continue
                self.successful_runs+=1
                lowest_rmsds.readline()
                rmsds = lowest_rmsds.readline().split(';')
                outfile.write(';'.join([case, rmsds[0], rmsds[1], rmsds[2]])+'\n')
                benchmark_results['lowest_10k'].append(float(rmsds[0]))
                benchmark_results['lowest_1k'].append(float(rmsds[1]))
                benchmark_results['lowest_10'].append(float(rmsds[2]))

        self.benchmark_results = benchmark_results

    def get_statistics(self):
        treshold = {'high quaility': 3.0, 'medium quality': 5.5}
        self.statistics = {}
        for type in ['10k', '1k', '10']:
            hq = len([rms for rms in self.benchmark_results['lowest_'+type] if rms < treshold['high quaility']])
            mq = len([rms for rms in self.benchmark_results['lowest_'+type] if treshold['high quaility'] <= rms <=
                                  treshold['medium quality']])
            lq = len([rms for rms in self.benchmark_results['lowest_'+type] if rms > treshold['medium quality']])
            self.statistics[type] = (hq, mq, lq)
        print self.statistics

    def print_summary(self):
        for (key, val) in self.statistics.items():
            message = 'Statistics for {}:\n'.format(key)+'High quality: {}\nMedium quality {}\nLow quality {}\n'.format(*[element/self.successful_runs for element in val])
            print message

    def sort_outfiles(self):
        self.plots_dir = '{}/plots'.format(self.done_benchdir)
        try:
            mkdir(self.plots_dir)
        except OSError:
            pass
        for dir in ['E_RMSD', 'RMSD_frame', 'RMSF', 'contact_maps', 'medoid_rmsds']:
            try:
                mkdir(self.plots_dir+'/'+dir)
            except OSError:
                pass

        for case in self.cases:
            for rzecz in ['E_RMSD', 'RMSD_frame', 'RMSF']:
                original_paths = glob('{}/run/{}/plots/{}*'.format(self.done_benchdir, case, rzecz))
                #print original_paths
                new_paths = ['{}/plots/{}/{}_'.format(self.done_benchdir, rzecz, case)+path.split('/')[-1] for path in original_paths]
                #print(new_paths)
                for org, nw in zip(original_paths, new_paths):
                    copyfile(org, nw)
                maps_original_paths = glob('{}/run/{}/contact_maps/*'.format(self.done_benchdir, case))
                maps_new_paths = ['{}/plots/{}/{}_'.format(self.done_benchdir, 'contact_maps', case)+path.split('/')[-1] for path in maps_original_paths]
                for org, nw in zip(maps_original_paths, maps_new_paths):
                    copyfile(org,nw)
                rmsd_medoids_original_paths = glob('{}/run/{}/output_data/medoids_rmsd*'.format(self.done_benchdir, case))
                #print rmsd_medoids_original_paths
                rmsd_medoids_new_paths = ['{}/plots/{}/{}_'.format(self.done_benchdir, 'medoid_rmsds', case)+path.split('/')[-1] for path in rmsd_medoids_original_paths]
                for org, nw in zip(rmsd_medoids_original_paths, rmsd_medoids_new_paths):
                    copyfile(org, nw)


#br = BenchmarkRunner(benchmark_file='./benchmark_data/2.txt')
#br = BenchmarkRunner(benchmark_file='./benchmark_data/MB_bench__bound_1.txt')
#br = BenchmarkRunner(benchmark_file='./kihara_peptides/kihara_bound.txt', runtype='frompdb')
#br.run_benchmark(test=True)
# # print br.benchdir
#runs = ['benchrun_Wed_Aug_30_17:15:59_2017', 'benchrun_Wed_Aug_30_17:18:42_2017', 'benchrun_Wed_Aug_30_17:19:06_2017']
#for run in runs:
#    ba = BenchmarkAnalyser('/Users/maciek/Desktop/cabsDock/'+run)
#    ba.sort_outfiles()
# ba.read_rmsds()
# ba.get_statistics()
# ba.print_summary()
# with open('benchmarking_logfile_{}.txt'.format(time.strftime('%x').replace('/','')), 'a+b') as log:
#     print("Bound benchmarks:")
#     log.write("Bound benchmarks:\n")
#     for i in xrange(3):
#         print('... Run {}'.format(i))
#         br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_bound_cases.txt', runtype='bound')
#         command = br.run_benchmark(test=True)
#         log.write(command)
#         log.write(';\n')
#     print("Unbound benchmarks:")
#     log.write("Unbound benchmarks:\n")
#     for i in xrange(3):
#         br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_unbound_cases.txt', runtype='unbound')
#         command = br.run_benchmark(test=True)
#         log.write(command)
#         log.write(';\n')
#ba = BenchmarkAnalyser('/Users/maciek/PyProjects/cabsdock/cabsDock/benchrun_Sat_Sep__9_00:14:45_2017')
#ba.get_pearsons()



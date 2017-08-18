from glob import glob
from os import mkdir, getcwd
from shutil import copyfile

from cabsDock.pbsgen import PbsGenerator
import time
from subprocess import call

class CommandGenerator(object):
    option_dictionary = {
        '-r': ('2gb1', 'aaaa'),
        '-p': None,
        '-e': None,
        '--excluding': None,
        '-f': None,
        '-R': None,
        '-P': None,
        '--separation': None,
        '--insertion-clash': None,
        '--insertion-attempts': None,
        '--ca-rest-add': None,
        '--sc-rest-add': None,
        '--ca-rest-weight': None,
        '-sc-rest-weight': None,
        '--ca-rest-file': None,
        '--sc-rest-file': None,
        '--mc-annealing': None,
        '--mc-cycles': None,
        '--mc-steps': None,
        '--replicas': None,
        '--replicas-dtemp': None,
        '--temperature': None,
        '-s': None,
        '--no-aa-rebuild': None,
        '--modeller-iterations': None,
        '--reference-pdb': None,
        '--clustering-iterations': None,
        '--filtering-number': None,
        '--clustering-medoids': None,
        '--load-cabs-files': None,
        '--contact-maps': None,
        '--align': None,
        '--reference-alignment': None,
        '--output-models': None,
        '--output-clusters': None,
        '--output-trajectories': None,
        '-c': None,
        '--image-file-format': None,
        '--work-dir': None,
        '--dssp-command': None,
        '--stride-command': None,
        '--fortran-command': None,
        '--save-config-file': None,
        '--save-cabs-files': None,
        '-V': None,

    }
    def __init__(self):
        pass


class BenchmarkRunner(object):
    def __init__(self, benchmark_file, options={'--image-file-format':'png'}, name=''):
        self.benchmark_file = benchmark_file
        self.options = options
        self.name = name

    def setup(self):
        if self.name == '':
            benchdir = getcwd()+'/benchrun_'
        else:
            benchdir = self.name
        benchdir+=time.strftime("%c").replace(' ','_')+''
        self.benchdir = benchdir
        try:
            mkdir(benchdir)
        except OSError:
            pass
        self.pbsgen = PbsGenerator(benchmark_list=self.benchmark_file, nonstandard_options_dict=self.options,
                                   rundir=self.benchdir+'/run')
        self.pbsgen.pbs_script(scriptdir=benchdir+'/pbs')

    def save_log(self):
        with open(self.benchdir+'/logfile', 'w') as logfile:
            logfile.write('#Run started: '+time.strftime("%c"))
            print self.options
            options_as_str = '\n'.join([str(key)+' = '+str(val) for (key,val) in self.options.items()])
            print options_as_str
            logfile.write('\n#Used options:\n'+options_as_str)
            cases_as_str = '\n'.join([case.work_dir for case in self.pbsgen.cases])
            logfile.write('\n#Cases:\n'+cases_as_str)

    def run_benchmark(self, test=False, qsub_options={'-l walltime':'60:00:00'}):
        self.setup()
        self.save_log()
        options_as_str = ' '.join([str(key)+'='+str(val) for (key, val) in qsub_options.items()])
        command = 'for f in {}/pbs/*.pbs; do qsub {} $f; done'.format(self.benchdir, options_as_str)
        if test:
            print command
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
                print line
                thisline = line.split()
                if 'Cases' in thisline[0]:
                    writecase = True
                elif writecase and '#' in thisline[0]:
                    writecase = False
                elif writecase:
                    self.cases.append(thisline[0])


    def read_rmsds(self):
        benchmark_results = {
            'lowest_10k': [],
            'lowest_1k': [],
            'lowest_10': [],
        }
        self.successful_runs = 0.
        with open(self.done_benchdir+'/lowest_rmsds_benchmark.txt') as outfile:
            for case in self.cases:
                lowest_rmsd_files = glob(self.done_benchdir + '/run/' + case + '/output_data/lowest_rmsds*')
                try:
                    lowest_rmsds = open(lowest_rmsd_files[0])
                except:
                    print('invalid rmsd file for case {}'.format(case))
                    continue
                self.successful_runs+=1
                lowest_rmsds.readline()
                rmsds = lowest_rmsds.readline().split(';')
                outfile.write(';'.join([case, rmsds[0], rmsds[1], rmsds[2]]))
                benchmark_results['lowest_10k'].append(rmsds[0])
                benchmark_results['lowest_1k'].append(rmsds[1])
                benchmark_results['lowest_10'].append(rmsds[2])

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

    def sort_pictures(self):
        self.plots_dir = '{}/plots'.format(self.done_benchdir)
        mkdir(self.plots_dir)
        for dir in ['E_RMSD', 'RMSD_frame', 'RMSF']:
            mkdir(self.plots_dir+'/'+dir)

        for case in self.cases:
            for rzecz in ['E_RMSD', 'RMSD_frame', 'RMSF']:
                copyfile(glob('{}/run/{}/plots/*.svg'.format(self.done_benchdir, case)), '{}/plots/{}'.format(self.done_benchdir, rzecz))



br = BenchmarkRunner(benchmark_file='./benchmark_data/2.txt')
#br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_cases.txt')
br.run_benchmark(test=True)
# # print br.benchdir
# ba = BenchmarkAnalyser('./benchbench')
# ba.read_rmsds()
# ba.get_statistics()
# ba.print_summary()

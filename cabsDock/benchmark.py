from glob import glob
from os import mkdir
from shutil import copyfile

from cabsDock.pbsgen import PbsGenerator
import time
from subprocess import call

class BenchmarkRunner(object):
    def __init__(self, benchmark_file, options={'--image-file-format':'png', '--work-dir':'./'}, name=''):
        self.benchmark_file = benchmark_file
        self.options = options
        self.name = name

    def setup(self):
        if self.name == '':
            benchdir = './benchrun_'
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
        options_as_str = ' '.join([str(key)+' '+str(val) for (key, val) in qsub_options.items()])
        command = 'for f in *.pbs; do qsub {} $f; done'.format(options_as_str)
        if test:
            print command
        else:
            call(command)

class BenchmarkAnalyser(object):
    def __init__(self, benchdict):
        self.benchdict = benchdict
        self.cases = []
        self.read_log()


    def read_log(self):
        benchlog = glob(self.benchdict+'/logfile')
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
        print self.cases
        benchmark_results = {
            'lowest_10k': [],
            'lowest_1k': [],
            'lowest_10': [],
        }
        self.successful_runs = 0.
        for case in self.cases:
            lowest_rmsd_files = glob(self.benchdict + '/' + case + '/output_data/lowest_rmsds*')
            print lowest_rmsd_files
            try:
                lowest_rmsds = open(lowest_rmsd_files[0])
            except:
                print('invalid rmsd file for case {}'.format(case))
                continue
            self.successful_runs+=1
            lowest_rmsds.readline()
            rmsds = lowest_rmsds.readline().split(';')
            benchmark_results['lowest_10k'].append(rmsds[0])
            benchmark_results['lowest_1k'].append(rmsds[1])
            benchmark_results['lowest_10'].append(rmsds[2])
        print benchmark_results
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
        mkdir('./plots')
        for dir in ['E_RMSD', 'RMSD_frame', 'RMSF']:
            mkdir('./plots'+dir)

        for case in self.cases:
            for rzecz in ['E_RMSD', 'RMSD_frame', 'RMSF']:
                copyfile(glob('./'+case+'/plots/'+'*.csv'), './plots/'+rzecz)



# br = BenchmarkRunner(benchmark_file='./benchmark_data/benchmark_cases.txt')
# br.run_benchmark(test=True)
# # print br.benchdir
# ba = BenchmarkAnalyser('./benchbench')
# ba.read_rmsds()
# ba.get_statistics()
# ba.print_summary()

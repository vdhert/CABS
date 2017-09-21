from os import mkdir
from os.path import exists, isdir


class PbsGenerator(object):
    def __init__(self, benchmark_list='', rundir='', nonstandard_options_dict={}, runtype='standard'):
        self.rundir = rundir
        self.benchmark_list = benchmark_list
        self.options_dict = nonstandard_options_dict
        self.cases = []
        f = open(self.benchmark_list, 'r')
        lines = f.readlines()

        for line in lines:
            row = line.split()
            if row[0] =='-':
                continue
            if runtype=='bound':
                #print row
                receptor = str(row[0]) + ':' + str(row[1])
                ligand = str(row[2] + ':' + str(row[3]))
                reference_pdb = str(row[0]) + ':' + str(row[1]) + str(row[4])
                sc_rests = []
                if len(row) > 5:
                    sc_rests = (str(row[5]), str(row[6]))
            elif runtype=='unbound':
                receptor = str(row[0]) + ':' + str(row[1])
                ligand = str(row[2] + ':' + str(row[3]))
                reference_pdb = str(row[4]) + ':' + str(row[5])+str(row[6])
                sc_rests = []
                #print row
                if len(row) > 7:
                    sc_rests = (str(row[7]), str(row[8]))
            elif runtype=='frompdb':
                receptor = str(row[0]) + ':' + str(row[1])
                ligand = str(row[0]) + ':' + str(row[2])
                reference_pdb = str(row[0]) + ':' + str(row[1]) + str(row[2])
                sc_rests = []
                #print row
                if len(row) > 3:
                    sc_rests = (str(row[3]), str(row[4]))

            elif runtype=='flex':
                receptor = str(row[0])
                ligand = None
                reference_pdb = None
                sc_rests = None

            self.cases.append(
                Case(
                    receptor = receptor,
                    ligand = ligand,
                    reference_pdb = reference_pdb,
                    sc_rests=sc_rests
                    )
                )

        if exists(self.rundir):
            if not isdir(self.rundir):
                raise Exception()
        else:
            mkdir(self.rundir)

        self.standard_err_out = '#PBS -o {}/out\n#PBS -e {}/err\n'
        self.standard_header = '#!/bin/bash\n'
        self.standard_cd = 'cd {}\n'.format(self.rundir)


    def pbs_script(self, scriptdir):
        if exists(scriptdir):
            if not isdir(scriptdir):
                raise Exception('File %s already exists and is not a directory' % scriptdir)
        else:
            mkdir(scriptdir)
        additional_options = ' '.join(['{} {}'.format(key, value) for (key, value) in self.options_dict.items()])
        for case in self.cases:
            name = case.receptor.split(':')[0]
            #print name
            with open(scriptdir+'/{}.pbs'.format(name), 'w') as scriptfile:
                scriptfile.write(self.standard_header)
                scriptfile.write(
                    self.standard_err_out.format(
                        self.rundir + '/{}'.format(name),
                        self.rundir + '/{}'.format(name)
                        )
                    )
                scriptfile.write(self.standard_cd)
                scriptfile.write(case.run_command()+' '+additional_options+'\n')

    def py_script(self, scriptdir):
        try:
            mkdir(scriptdir)
        except OSError:
            pass
        for case in self.cases:
            name = case.receptor.split(':')[0]
            # print name
            with open(scriptdir + '/{}.py'.format(name), 'w') as pyfile:
                pyfile.write('from cabsDock.job import FlexTask\n')
                command = 'flextask = FlexTask('
                command += 'structure=\'{}\', temperature=(1.4, 1.4), receptor_restraints=(\'ss2\', 3, 3.8, 8.0),'.format(name)
                command += ' mc_annealing=20, mc_cycles=50, mc_steps=50, work_dir=\'{}/flextask_{}\')'.format(self.rundir, name)
                command += '\nflextask.run()'
                print command
                pyfile.write(command)
            with open(scriptdir+'/{}.pbs'.format(name), 'w') as pbsfile:
                pbsfile.write(self.standard_header)
                # pbsfile.write(
                #     self.standard_err_out.format(
                #         self.rundir + '/{}'.format(name),
                #         self.rundir + '/{}'.format(name)
                #     )
                # )
                pbsfile.write(self.standard_cd)
                pbsfile.write('python ../pbs/{}.py'.format(name))


class Case(object):
    def __init__(
            self,
            receptor=None,
            ligand=None,
            reference_pdb=None,
            sc_rests=[],
            ):
        super(Case, self).__init__()
        self.receptor = receptor
        self.ligand = ligand
        self.work_dir = str(receptor.split(':')[0])
        self.reference_pdb = reference_pdb
        self.sc_rests = sc_rests

    def __str__(self):
        return ' '.join([str(self.receptor), str(self.ligand), str(self.work_dir), str(self.reference_pdb)])

    def run_command(self):
        command = 'cabsDock -r {} -p {} --work-dir {} --reference-pdb {}'.format(self.receptor, self.ligand, self.work_dir, self.reference_pdb)
        if self.sc_rests:
            command += ' --sc-rest-add {}'.format(self.sc_rests[0]+' '+self.sc_rests[1]+' 5.0 1.0')
        return command
#usage
# pbsgntr = PbsGenerator(benchmark_list='./benchmark_data/cabsflex_onechain.txt', runtype='flex', rundir='.')
# pbsgntr.py_script('pys')
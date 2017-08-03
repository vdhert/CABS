from os import mkdir
from os.path import exists, isdir


class PbsGenerator(object):
    def __init__(self, benchmark_list='', nonstandard_options_dict={}):
        self.benchmark_list = benchmark_list
        self.options_dict = nonstandard_options_dict
        self.cases = []
        f = open(self.benchmark_list, 'r')
        lines = f.readlines()


        for line in lines:
            row = line.split()
            self.cases.append(Case(
                    receptor = str(row[0])+':'+str(row[1]),
                    ligand = str(row[3]+':'+str(row[4])),
                    reference_pdb = str(row[0])

                )
            )

        self.standard_header='#!/bin/bash\ncd ~/cabsDock/benchmark/\n'

    def pbs_script(self, scriptdir):
        if exists(scriptdir):
            if not isdir(scriptdir):
                raise Exception('File %s already exists and is not a directory' % work_dir)
        else:
            mkdir(scriptdir)
        additional_options = ' '.join(['{} {}'.format(key, value) for (key, value) in self.options_dict.items()])
        for case in self.cases:
            with open(scriptdir+'/{}.pbs'.format(case.receptor.split(':')[0]), 'w') as scriptfile:
                scriptfile.write(self.standard_header)
                scriptfile.write(case.run_command()+' '+additional_options+'\n')

class Case(object):
    def __init__(
            self,
            receptor=None,
            ligand=None,
            reference_pdb=None,
            ):
        super(Case, self).__init__()
        self.receptor = receptor
        self.ligand = ligand
        self.work_dir = str(receptor.split(':')[0])
        self.reference_pdb = reference_pdb

    def __str__(self):
        return str(self.receptor)+' '+str(self.ligand)+' '+str(self.work_dir)+' '+str(self.reference_pdb)

    def run_command(self):
        return 'cabsDock -r {} -p {} --work-dir {} --reference-pdb {}'.format(self.receptor, self.ligand, self.work_dir, self.reference_pdb)

#usage
# pbsgntr = PbsGenerator(benchmark_list='./benchmark_data/benchmark_cases.txt')
# pbsgntr.pbs_script('pbs')
import csv
import argparse
import re
from pkg_resources import require

__version__ = require('CABS')[0].version


def mk_parser(def_prot_rest, def_temp, def_replicas):
    parser = argparse.ArgumentParser()
    parser.add_argument_group('BASIC OPTIONS')
    parser.add_argument(    '-P',
                            '--protein',
                            )

    parser.add_argument_group('PROTEIN OPTIONS')
    parser.add_argument(    "-f",
                            "--protein-flexibility",
                            metavar="FLEXBILITY",
                            default=1.0,
                            type=float,
                            help="""""")

    parser.add_argument(    "-R",
                            "--protein-restraints",
                            metavar="MODE GAP MIN MAX",
                            nargs=4,
                            default=def_prot_rest,
                            help="""""")

    parser.add_argument_group('RESTRAINTS OPTIONS', "Restraints options allows to set up distance restraints, either between C-alpha atoms (CA) or Side Chains (SC), where SCs are geometric centers of side chains atoms (as defined in CABS coarse-grained model).")
    parser.add_argument(    "--ca-rest-add",
                            metavar="RESI RESJ DIST WEIGHT",
                            action='append',
                            help="""Adds distance restraint between C-alpha (CA) atom in residue RESI and C-alpha atom in residue RESJ").

                                DIST is a distance between these atoms and WEIGHT is restraint's weight from [0, 1].

                                In order to add restraints between the peptide and the receptor, or between two peptides, use PEP1, PEP2, ... as chain identifiers of the peptides (even when peptide is read from a pdb file its chain identifier is ignored).

                                i.e. ""123:A 5:PEP1 8.7 1.0"" adds restraint between the C-alpha atom of the residue number 123 in the chain A of the receptor and the C-alpha atom of the 5th residue of the peptide.

                                If you add only one peptide both 'PEP' and 'PEP1' is a valid chain identifier.

                                If you add multiple peptides they will be ordered as follows: [1] from config file added by the 'peptide' option [2] from config file added by the 'add-peptide' option [3] from command line added by the '-p, --peptide' option [4] from command line added by the '--add-peptide' option.

                                Peptides added by the same method preserve the order by which they appear in the config file, or on the command line.

                                Can be used multiple times to add multiple restraints."
                        """)

    parser.add_argument(    "--sc-rest-add",
                            metavar="RESI RESJ DIST WEIGHT",
                            action='append',
                            help="""Adds distance restraint between SC pseudoatom in residue RESI and SC pseudoatom in residue RESJ; DIST is a distance between these pseudoatoms (geometric centers of side chain atoms) and WEIGHT is restraint's weight from [0, 1]. Can be used multiple times to add multiple restraints.""")

    parser.add_argument(    "--ca-rest-weight",
                            metavar="WEIGHT",
                            action='append',
                            default=1.0,
                            type=float,
                            help="""Set global weight for all C-alpha restraints (including automatically generated restraints for the receptor).""")

    parser.add_argument(    "--sc-rest-weight",
                            metavar="WEIGHT",
                            action='append',
                            default=1.0,
                            type=float,
                            help="""Set global weight for all SC restraints.""")

    parser.add_argument(    "--ca-rest-file",
                            metavar="FILE",
                            action='append',
                            type=open,
                            help="""Read C-alpha restraints from file (use multiple times to add multiple files).""")

    parser.add_argument(    "--sc-rest-file",
                            metavar="FILE",
                            action='append',
                            type=open,
                            help="""Read SC restraints from file (use multiple times to add multiple files).""")

    parser.add_argument_group('SIMULATION OPTIONS', description="""Simulation options allow to modify different parameters of Replica Exchange Monte Carlo simulation procedure""")

    parser.add_argument(    "--mc-annealing",
                            metavar="NUM",
                            type=int,
                            default=20,
                            help="""Sets number of Monte Carlo temperature annealing cycles to NUM (NUM > 0, default value = 20, changing default value is recommended only for advanced users).""")

    parser.add_argument(    "--mc-cycles",
                            metavar="NUM",
                            type=int,
                            default=50,
                            help="""Sets number of Monte Carlo cycles to NUM (NUM>0, default value = 50). Total number of snapshots generated for each replica/trajectory = [mc-annealing] x [mc-cycles], default: 20x50=1000.""")

    parser.add_argument(    "--mc-steps",
                            metavar="NUM",
                            type=int,
                            default=50,
                            help="""Sets number of Monte Carlo cycles between trajectory frames to NUM (NUM > 0, default value = 50). NUM = 1 means that every generated conformation will occur in trajectory. This option enables to increase the simulation length (between printed snapshots) and doesn't impact the number of snapshots in trajectories (see also --mc-steps description).""")

    parser.add_argument(    "--replicas",
                            metavar="NUM",
                            type=int,
                            default=def_replicas,
                            help="""Sets number of replicas to be used in Replica Exchange Monte Carlo (NUM > 0, default value = 10, changing default value is recommended only for advanced users).""")

    parser.add_argument(    "--replicas-dtemp",
                            metavar="DELTA",
                            type=float,
                            default=0.5,
                            help="""Sets temperature increment between replicas (DELTA > 0, default value = 0.5, changing default value is recommended only for advanced users).""")

    parser.add_argument(    "--temperature",
                            metavar="TINIT TFINAL",
                            nargs=2,
                            type=float,
                            default=def_temp,
                            help="""Sets temperature range for simulated annealing TINIT - initial temperature, TFINAL - final temperature (default values TINIT=2.0, TFINAL=1.0, changing default value is recommended only for advanced users). CABSdock uses temperature units, which do not correspond straightforwardly to real temperatures. Temperature around 1.0 roughly corresponds to nearly frozen conformation, folding temperature of small proteins in the CABS model is usually around 2.0.""")

    parser.add_argument(    "-s",
                            "--random-seed",
                            metavar="NUM",
                            type=int,
                            help="""Sets seed for random number generator.""")

    parser.add_argument_group('ALL-ATOM RECONSTRUCTION OPTIONS', description="""All-atom reconstruction options allow to set up details of all-atom reconstruction and refinement procedure, which is performed using MODELLER package""")

    parser.add_argument(    "--no-aa-rebuild",
                            action='store_true',
                            help="""Skip the all-atom reconstruction. Saves CA models for the top scored binding poses instead.""")

    parser.add_argument(    "--modeller-iterations",
                            metavar='NUM',
                            default=3,
                            type=int,
                            help="""Set number of iterations for reconstruction procedure in MODELLER package (default: 3). Bigger numbers may result in more accurate models, but reconstruction takes longer.""")

    parser.add_argument_group('ANALYSIS OPTIONS', description='''Analysis options allow to perform comparison analyses (to provided reference complex structure) and for repeated scoring and analysis of CABSdock trajectories.''')

    parser.add_argument(    "--reference-pdb",
                            metavar='PDB',
                            default=3,
                            type=int,
                            help="""Load reference complex structure. This option allows for comparison with the reference complex structure and triggers additional analysis features

                            PDB must be either:

                            * [pdb code]:[receptor chains]:[peptide(s) chain(s)]
                            * [pdb file]:[receptor chains]:[peptide(s) chain(s)]

                            i.e 1abc:AB:C, 1abc:AB:CD, myfile.pdb:AB:C, myfile.pdb.gz:AB:CDE""")

    parser.add_argument(    "--clustering-iterations",
                            metavar='NUM',
                            default=100,
                            type=int,
                            help="""Sets number of iterations of the clustering k-medoids algorithm.""")

    parser.add_argument(    "--filtering-number",
                            metavar='NUM',
                            default=1000,
                            type=int,
                            help="""Sets number of low-energy models from trajectories to be clustered (default 1000).""")

    parser.add_argument(    "--filtering-fromeach",
                            action='store_true',
                            help="""Picks (filtering-number/replicas) models from each replica. If False picks (filtering-number) low-energy models from the whole trajectory.""")

    parser.add_argument(    "--clustering-medoids",
                            metavar='NUM',
                            default=10,
                            type=int,
                            help="""Sets number of medoids in k-medoids clustering algorithm.""")

    parser.add_argument(    "--load-cabs-files",
                            metavar='FILE',
                            type=int,
                            help="""Loads CABSdock simulation files and allows for repeated scoring and analysis of CABSdock trajectories (with new settings, for example using a reference complex structure and --reference option).""")

    parser.add_argument(    "--contact-maps",
                            action='store_true',
                            help="""Stores contact maps matrix plots and histograms of contact frequencies.""")

    parser.add_argument(    "--align",
                            metavar='METHOD',
                            default='SW',
                            help="""Method to be used to align terget with reference pdb. Available options are: SW -- Smith-Waterman (default), blastp -- protein BLAST (requires NCBI+ package installed), trivial -- simple sequential alignment, used only in case of one-chain input and reference of the same length.""")

    parser.add_argument(    "--alignment-options",
                            metavar='KEY1=VAL1 KEY2=VAL2...',
                            action='append',
                            default=[],
                            type=lambda x: x.split('='),
                            help="""Path to alignment with reference structure. If set 'align' option is ignored.""")

    parser.add_argument(    "--alignment-peptide-options",
                            metavar='KEY1=VAL1 KEY2=VAL2...',
                            action='append',
                            default=[],
                            type=lambda x: x.split('='),
                            help="""Path to alignment with reference structure. If set 'align' option is ignored.""")

    parser.add_argument(    "--cc-threshold",
                            metavar='NUM',
                            default=6.5,
                            type=float,
                            help="""Cmap contact criterion for side chain (SC) distance threshold.""")

    parser.add_argument_group('OUTPUT OPTIONS', description='Output options')

    parser.add_argument(    "--output-models",
                            metavar='NUM',
                            default=10,
                            type=int,
                            help="""Set number of final models to be generated.""")

    parser.add_argument(    "--output-clusters",
                            action='store_true',
                            help="""Save pdb files with clusters.""")

    parser.add_argument(    "--output-trajectories",
                            action='store_true',
                            help="""Save pdb files with trajectories.""")

    parser.add_argument_group('MISCELLANEOUS OPTIONS', description='Miscellaneous options')

    parser.add_argument(    "-c",
                            "--config",
                            metavar='CONFIG',
                            type=open,
                            help="""Read options from configuration file CONFIG.""")

    parser.add_argument(    "--image-file-format",
                            metavar='FMT',
                            default='svg',
                            help="""Produces all the image files in given format.""")


    parser.add_argument(    "--work-dir",
                            metavar='DIR',
                            default='.',
                            help="""""")

    parser.add_argument(    "--dssp-command",
                            metavar='PATH',
                            default='mkdssp',
                            help="""""")

    parser.add_argument(    "--fortran-command",
                            metavar='PATH',
                            help="""""")

    parser.add_argument(    "--save-config-file",
                            metavar='PATH',
                            default='svg',
                            help="""""")

    parser.add_argument(    "--save-cabs-files",
                            action='store_true',
                            help="""""")

    parser.add_argument(    "-V",
                            "--verbose",
                            metavar='VERBOSITY',
                            default=1,
                            type=int,
                            help="""""")

    return parser

def mk_flex_parser():
    parser = mk_parser(def_prot_rest=('ss2', 3, 3.8, 8.0), def_temp=(1.4, 1.4), def_replicas=1)
    return parser

def mk_dock_parser():
    parser = mk_parser(def_prot_rest=('all', 5, 5.0, 15.0), def_temp=(2.0, 1.0), def_replicas=10)
    return parser

class ConfigFileParser:

    OPTIONRE = re.compile(
        r'(?P<option>[^:=]*)'
        r'[:=]'
        r'(?P<value>.*)$'
    )

    def __init__(self, filename):
        self.args = []
        with open(filename) as f:
            for line in f:
                if line == '' or line[0] in ';#\n':
                    continue
                match = self.OPTIONRE.match(line)
                if match:
                    option, value = match.groups()
                    self.args.append('--' + option.strip())
                    self.args.extend(value.split('#')[0].split(';')[0].split())
                else:
                    self.args.extend(line.split('#')[0].split(';')[0].split())


if __name__ == '__main__':
    pass

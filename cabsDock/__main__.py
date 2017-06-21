import argparse

intro="""
cabsDock - protein peptide docking tool
"""

goodbye="""
visit biocomp.chem.uw.edu.pl/CABSDock
"""


def run_job():
    parser = argparse.ArgumentParser(
        prog='cabsDock',
        description=intro,
        epilog=goodbye,
        formatter_class=argparse.RawTextHelpFormatter,
        version=0.99
    )

    parser.add_argument(
        'receptor',
        metavar = 'RECEPTOR',
        help="""Load receptor structure.
%(metavar)s must be a pdb file or a pdb_code
(optionally with chain id(s) i.e. 1rjk:AC)
"""
    )
    
    parser.add_argument(
        '-p', '--add-peptide',
        nargs='+',
        action='append',
        dest='ligand',
        metavar=('PEPTIDE', 'CONFORMATION LOCATION'),
        help="""
Add peptide to the complex.
This option can be used multiple times to add multiple peptides.

SEQUENCE must be either:
    1. SEQUENCE in one letter code (optionally with SECONDARY STRUCTURE
       C - coil, H - helix, E - sheet, T - turn, i.e. LAHCIM:CHHTEC)
    2. pdb file (may be gzipped)
    3. pdb code (optionally with chain id i.e 1rjk:C)

CONFORMATION sets initial conformation of the peptide. Must be either:
    1. random - conformation is randomized
    2. keep - preserve conformation from file.
    This has no effect if PEPTIDE=SEQUENCE
    
LOCATION sets initial location for the peptide. Must be either:
    1. random - peptide is placed in a random location
       at 20A(default value, can be changed)
       from the receptor's surface
    2. keep - preserve location from file.
    This has no effect if PEPTIDE=SEQUENCE
    3. patch - where patch is a list of receptor residues
       joined with + specifying where to put the peptide (i.e 123:A+125:A+17:B)
    WARNING: residues in patch should be on the surface
    of the receptor and close to each other.
"""
    )
    parser.add_argument(
        '-d', '--dir',
        dest='work_dir',
        help='Set working directory to WORK_DIR. Default is current dir.'
    )
    parser.add_argument(
        '-c', '--config',
        help='read options from CONFIG file'
    )
    parser.add_argument(
        '-n', '--mc_cycles',
        type=int,
        help='number of MC cycles, controls number of models in a trajectory'
    )
    parser.add_argument(
        '-m', '--mc_steps',
        type=int,
        help='number of MC steps between trajectory snapshots'
    )
    parser.add_argument(
        '-ti', '--initial_temperature',
        type=float,
        dest='t_init',
        help='initial temperature'
    )
    parser.add_argument(
        '-tf', '--final_temperature',
        type=float,
        dest='t_final',
        help='final temperature'
    )
    parser.add_argument(
        '-nr', '--replicas',
        type=int,
        help='number of replicas'
    )
    parser.add_argument(
        '-dssp', '--dssp_command',
        help='command used to run DSSP'
    )

    args = parser.parse_args()
    print args
    #job_args = {k: v for k, v in vars(args).items() if v}
    #from job import Job
    #j = Job(**job_args)
    #j.run_job()

if __name__ == '__main__':
    run_job()

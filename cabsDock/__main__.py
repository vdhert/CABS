import argparse

intro="""
cabsDock - protein peptide docking tool
"""

goodbye="""
visit biocomp.chem.uw.edu.pl/CABSDock
"""

add_peptide_help="""
Add peptide to the complex.
This option can be used multiple times to add multiple peptides.

PEPTIDE must be either:

    1. SEQUENCE in one-letter code
       optionally followed by :SECONDARY_STRUCTURE (i.e. LAHCIM:CHHTEC)
       where C - coil, H - helix, E - sheet, T - turn

    2. pdb file (may be gzipped)

    3. pdb code (optionally with chain id i.e 1abc:D)

CONFORMATION sets initial conformation of the peptide. Must be either:

    1. random - conformation is randomized [default]

    2. keep - preserve conformation from file
       This has no effect if PEPTIDE=SEQUENCE.
    
LOCATION sets initial location for the peptide. Must be either:

    1. random - peptide is placed in a random location on the surface
       of a sphere centered at the receptor's geometrical center
       at distance defined by the '-s, --separation' option
       from the surface of the receptor.
       
    2. keep - preserve location from file.
       This has no effect if PEPTIDE=SEQUENCE

    3. patch - list of receptor's residues (i.e 123:A+125:A+17:B)
       Peptide will be placed above the geometrical center of listed
       residues at distance defined by the '-s, --separation' option
       from the surface of the receptor.
       WARNING: residues listed in path should be on the surface
       of the receptor and close to each other.
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
        '-r', '--receptor',
        metavar = 'RECEPTOR',
        dest='receptor',
        required=True,
        help="""Load receptor structure.
%(metavar)s must be a pdb file or a pdb_code
(optionally with chain id(s) i.e. 1abc:DE)
"""
    )
    
    parser.add_argument(
        '-p', '--add-peptide',
        nargs='+',
        action='append',
        dest='ligand',
        metavar=('PEPTIDE', 'CONFORMATION LOCATION'),
        help=add_peptide_help
    )
    parser.add_argument(
        '-d', '--dir',
        dest='work_dir',
        metavar='DIR',
        help='Set working directory to WORK_DIR. Default is current dir.'
    )
    parser.add_argument(
        '-c', '--config',
        dest='config',
        metavar='CONFIG',
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

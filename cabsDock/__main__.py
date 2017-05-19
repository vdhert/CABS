import argparse


def run_job():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-r', '--receptor',
        help='Load receptor structure.\n'
             '[receptor] must be a pdb file or a pdb_code\n'
             '(optionally with chain id(s) i.e. 1rjk:AC)',
        required=True
    )
    parser.add_argument(
        '-p', '--add_peptide',
        nargs=3,
        action='append',
        dest='ligands',
        metavar=(
            'PEPTIDE',
            'CONFORMATION',
            'LOCATION'
        ),
        help='Add peptide to the complex.\n'
             'This option can be used multiple times to add multiple peptides.\n'
             'PEPTIDE must be either:\n'
             '1.\tSEQUENCE in one letter code\n'
             '\t(optionally with SECONDARY STRUCTURE\n'
             '\tC - coil, H - helix, E - sheet, T - turn\n'
             '\ti.e. LAHCIM:CHHTEC)\n'
             '2.\tpdb file\n'
             '3.\tpdb_code (optionally with chain id i.e 1rjk:C)\n'
             '\nCONFORMATION sets initial conformation of the peptide.\n'
             'Must be either:\n'
             '1.\trandom - conformation is randomized\n'
             '2.\tkeep - preserve conformation from file.\n'
             '\tThis has no effect if PEPTIDE=SEQUENCE\n'
             '\nLOCATION sets initial location for the peptide.\n'
             'Must be either:\n'
             '1.\trandom - peptide is placed in random location\n'
             '\tat 20A(default value, can be changed) from the receptor\'s surface\n'
             '2.\tkeep - preserve location from file.\n'
             '\tThis has no effect if PEPTIDE=SEQUENCE\n'
             '3.\tpatch - where patch is a list of receptor residues joined with \'+\'\n'
             '\tspecifying where to put the peptide (i.e 123:A+125:A+17:B)\n'
             '\tWARNING: residues in patch should be on the surface\n'
             '\tof the receptor and close to each other.\n'
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

    args = parser.parse_args()
    job_args = {k: v for k, v in vars(args).items() if v}
    from job import Job
    j = Job(**job_args)

if __name__ == '__main__':
    run_job()

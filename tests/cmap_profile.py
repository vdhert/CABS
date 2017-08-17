from cabsDock.job import Job
from argparse import ArgumentParser as AP


def run():
    j = Job(
            work_dir='./profile_res/',
            receptor='1klu:AB',
            ligand=[['GELIGTLNAAKVPAD:CCCEEEECCEECCCC', 'random', 'random']],
            mc_cycles=20,
            mc_steps=50,
            replicas=10,
            load_cabs_files=(cla.data_path + '1klu/TRAF', cla.data_path + '1klu/SEQ'),
            AA_rebuild=False
            )
    j.cabsdock()

if __name__ == '__main__':
    ap = AP()
    ap.add_argument('--data_path', default='./tests/data/')
    cla, args = ap.parse_known_args()
    run()

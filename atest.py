from cabsDock.job import Job
if __name__ == '__main__':
    j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=50,  mc_steps=1, replicas=10, dbg=True)

from cabsDock.job import Job
if __name__ == '__main__':
    j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=1,  mc_steps=1, replicas=1, dbg=True)
    j.run_job()

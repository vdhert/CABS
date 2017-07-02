from cabsDock.job import Job
if __name__ == '__main__':
    j = Job(receptor='1ce1:HL ', ligand=[['GTSSPSAD'],['AAAAAAA']], mc_cycles=1,  mc_steps=1, replicas=1, dbg=True)
    j.run_job()
    # j = Job(receptor='2gb1', ligand=[['MICHAL']], mc_cycles=2,  mc_steps=1, replicas=10, dbg=True)
    # j.run_job()

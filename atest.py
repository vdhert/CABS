from cabsDock.job import Job
if __name__ == '__main__':
    j = Job(
        receptor='1ce1:HL ',
        ligand=[['GTSSPSAD']],
        mc_cycles=5,
        mc_steps=5,
        replicas=10,
        native_pdb='1ce1',
        native_receptor_chain='HL',
        native_peptide_chain='P',
        dbg=True,
        work_dir='testdir_mcc5_mcs5_reps10')
    j.run_job()
    # j = Job(receptor='2gb1', ligand=[['MICHAL']], mc_cycles=2,  mc_steps=1, replicas=10, dbg=True)
    # j.run_job()

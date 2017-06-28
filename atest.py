from cabsDock.job import Job
if __name__ == '__main__':
    j = Job(receptor='1jbu:H', ligand = [['EEWEVLCWTWETCER']], mc_cycles=20, mc_steps=1, replicas=10, native_pdb='1jbu',
                               native_receptor_chain='H',
                               native_peptide_chain='X')
    print j.run_job()

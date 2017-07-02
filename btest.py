import pickle, numpy
from cabsDock.job import Job
from cabsDock.protein import ProteinComplex

try:
    from pymol import cmd
except: pass

with open('traj.pck') as f:
    tr = pickle.load(f)
with open('flti.pck') as f:
    flt = pickle.load(f)
with open('clst.pck') as f:
    c = pickle.load(f)

j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=2,  mc_steps=1, replicas=2, dbg=True)
j.initial_complex = ProteinComplex(j.config)
#~ tr.coordinates = numpy.array([tr.coordinates[6,:,:,:]])
j.mk_cmaps(tr, c, flt, 4.5)
import pickle, numpy
from cabsDock.job import Job
from cabsDock.protein import ProteinComplex
from cabsDock.utils import SCModeler
from cabsDock.cmap import ContactMapFactory

try:
    from pymol import cmd
except: pass

with open('traj.pck') as f:
    tr = pickle.load(f)
with open('flti.pck') as f:
    flt = pickle.load(f)
with open('clst.pck') as f:
    c = pickle.load(f)
with open("cplx.pck", ) as f:
    cplx = pickle.load(f)

tr.coordinates = tr.coordinates[:1, -100:, :68, :3]

scmodeler = SCModeler(cplx)
sc_traj_full = scmodeler.calculate_sc_traj(tr.coordinates)
print 'sc_traj'
targ_cmf = ContactMapFactory(cplx.receptor_chains, 'C', tr.template)
print 'cmf'
cmap0t = targ_cmf.mk_cmap(sc_traj_full, 20.)[0]
print 'cm'
cmap0t.save_histo('testhisto')
print 'histo'

#~ j = Job(receptor='2gb1', ligand=[['MICHAL'], ['LAHCIM']], mc_cycles=2,  mc_steps=1, replicas=2, dbg=True)

#~ j.initial_complex = ProteinComplex(j.config)
#~ tr.coordinates = numpy.array([tr.coordinates[6,:,:,:]])
#~ j.mk_cmaps(tr, c, flt, 4.5)
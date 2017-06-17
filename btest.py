import pickle
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap

from cabsDock.utils import SCModeler

try:
    from pymol import cmd
except: pass

with open('test_traj.pck') as f:
    trajectory = pickle.load(f)
with open('test_complex.pck') as f:
    initial_complex = pickle.load(f)
with open('test_clusters.pck') as f:
    C = pickle.load(f)

tra = trajectory.filter(1000)
lig = tra.select('chain B')
scm = SCModeler(trajectory.template)

for lig in initial_complex.ligand_chains:
    cmf = ContactMapFactory(initial_complex.receptor_chains, lig, trajectory.template, mth=scm.execute)      #where from comes rec ch name?
    cmaps = cmf.mk_cmap(trajectory.coordinates, 6.5)
    for n, cmap in enumerate(cmaps):
        with open('res_btest/cm_%s_%i.txt' % (lig, n), 'w') as f:
            cmap.save_txt(f)
        cmap.save_png('res_btest/cm_%s_%i' % (lig, n))
    #~ import pdb; pdb.set_trace()

#~ def fnc(v, n):
    #~ for j, i in enumerate(v): cmd.pseudoatom(n + str(j), pos=i)

#~ from string import ascii_letters

#~ fnc(map(tuple, thx), 'a')
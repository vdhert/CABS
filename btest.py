import pickle
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap

from cabsDock.utils import rebuildCb

with open('test_traj.pck') as f:
    trajectory = pickle.load(f)
with open('test_complex.pck') as f:
    initial_complex = pickle.load(f)
with open('test_clusters.pck') as f:
    C = pickle.load(f)

tra = trajectory.filter(1000)
lig = tra.select('chain B')

for lig in initial_complex.ligand_chains:
    cmf = ContactMapFactory(initial_complex.receptor_chains, lig, trajectory.template)      #where from comes rec ch name?
    rebuildCb(trajectory.coordinates[0, 0, cmf.inds1,], cmf.ats1)
    #~ cmaps = cmf.mk_cmap(trajectory.coordinates, 6.5)
    #~ for n, cmap in enumerate(cmaps):
        #~ with open('res_btest/cm_%s_%i.txt' % (lig, n), 'w') as f:
            #~ cmap.save_txt(f)
        #~ cmap.save_png('res_btest/cm_%s_%i' % (lig, n))
    #~ import pdb; pdb.set_trace()


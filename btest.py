import pickle
from cabsDock.cmap import ContactMapFactory
from cabsDock.cmap import ContactMap

with open('test_traj.pck') as f:
    trajectory = pickle.load(f)
with open('test_complex.pck') as f:
    initial_complex = pickle.load(f)
with open('test_clusters.pck') as f:
    C = pickle.load(f)

tra = trajectory.filter(1000)
lig = tra.select('chain B')

#parts to add to Job.__init__
#~ for lig in initial_complex.ligand_chains:
    #~ cmf = ContactMapFactory(initial_complex.receptor_chains, lig, trajectory.template)      #where from comes rec ch name?
    #~ cmaps = cmf.mk_cmap(trajectory.coordinates, 6.5)
    #~ for n, cmap in enumerate(cmaps):
        #~ cmap.save('btest_res/cm_%s_%i' % (lig, n))

import pdb; pdb.set_trace()

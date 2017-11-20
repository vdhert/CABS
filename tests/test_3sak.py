from CABS.pdblib import Pdb
from CABS import utils, logger
import numpy as np

logger.setup(4, remote=True)

reference = Pdb('3sak', selection='name CA and model 1').atoms
trajectory = Pdb('output_pdbs/replica.pdb').atoms

target = reference.to_numpy()

for _index, model in enumerate(trajectory.models(), 1):
    query = model.to_numpy()
    rmsd, rotation, target_shift, query_shift = utils.dynamic_kabsch(target, query)
    query = np.dot(query - query_shift, rotation) + target_shift
    model.from_numpy(query)
    print 'Model: %i Weighted rmsd: %8.3f Rmsd: %8.3f' % (_index, rmsd, utils.rmsd(target, query))

with open('super.pdb', 'wb') as f:
    f.write(trajectory.make_pdb())
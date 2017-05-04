from trajectory import Trajectory
from utilities import Utilities
from cluster import Clustering, Cluster
from cabsDock import pdb, atom
import numpy
import timeit
atoms_input = atom.Atoms(pdb.Pdb("Test_case/input.pdb"))
headers = numpy.load("Test_case/headers_stored.npy")

trajectories = [ Trajectory(atom.Atoms(pdb.Pdb("Test_case/trajectory_"+str(i+1)+".pdb")), "A", "B", 1, headers[i]) for i in xrange(10) ]
#trajectories = [ Trajectory(atom.Atoms(pdb.Pdb("Test_case/trajectory_"+str(i+1)+".pdb")), "A", "B", 1, headers[i]) for i in xrange(1) ]
input_structure = Trajectory(atoms_input.select("name CA"), "A", "B", 1, [ [0, 0, 0, 0, 0] ])

for trajectory in trajectories:
	trajectory = Utilities.fit_models(trajectory, input_structure)

list_of_sorted_models, filtered_atoms = Utilities.filter_models(trajectories = trajectories)
trajectory_ready_for_clustering = Trajectory(filtered_atoms, "A", "B", 1, [ [0, 0, 0, 0, 0] for i in range(len(list_of_sorted_models)) ])
clustering = Clustering(trajectory = trajectory_ready_for_clustering, input_model = input_structure)
medoids, clusters_list, labels = clustering.perform_clustering()

print(medoids)
print(clusters_list)

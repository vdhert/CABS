from cabsDock import pdb, atom
from utilities import Utilities
import numpy
from matplotlib import pyplot

class Cluster(object):
    """ Class for storing one cluster resulting from Clustering.perform_clustering() call. """
    def __init__(self, representative_structure = None, structures = [], feature = None, cluster_distance_matrix = None):
        super(Cluster, self).__init__()

        # representative structure of the cluster (i.e. medoid)
        self.representative_structure = representative_structure

        # list of structures in this cluster
        self.structures = structures

        # a ranking-related feature of this cluster (i.e. cluster density)
        self.feature = feature

        # a cluster distance matrix (submatrix of the distance matrix used to perform the clustering)
        self.cluster_distance_matrix = cluster_distance_matrix

    def _calculate_cluster_density(self, method = 'default'):
        # this method updates self.feature attribute AND returns its value

        # chosen measure of cluster volume (default: maximum distance between cluster members)
        cluster_volume_measure = { 'default' : max(self.cluster_distance_matrix)}

        # chosen measure of cluster mass (default: the number of structures in the cluster)
        cluster_mass_measure = { 'default' : len(self.structures)}

        self.feature = cluster_mass_measure[method] / cluster_volume_measure[method]
        return self.feature


class Clustering(object):
    """ Class for performing structural clustering of CABS-dock models and storing the results of clustering. """
    def __init__(self, trajectory = None, input_model = None):
        super(Clustering, self).__init__()
        # trajectory: Trajectory instance with models to cluster, input_model: Trajectory instance with one input model.
        self.trajectory = trajectory
        self.input_model = input_model

    def _calculateDistanceMatrix(self, structures_aligned = True): 
        """ An in-place method that updates self.distanceMatrix of a Clustering instance """
        if structures_aligned is not True:
            # list of all models aligned to the input structure
            list_of_models_to_cluster = Utilities.fit_models(self.trajectory, self.input_model)
        else:
            list_of_models_to_cluster = [ model for model in self.trajectory ] 
        # reading the ligand ids (Ogarnac wieksza liczbe ligandow)
        trajectory_ligand_id = self.trajectory.ligand_id[0]
        input_ligand_id = self.input_model.ligand_id

        #restricting the trajectories to CAs of ligands ( poprawic dla wiekszej liczby lancuchow receptora )
        self.trajectory.atoms = self.trajectory.atoms.select("chain "+str(trajectory_ligand_id)+" and name CA")
        self.input_model.atoms = self.input_model.atoms.select("chain "+str(input_ligand_id)+" and name CA")

        self.distanceMatrix = numpy.zeros(shape =  ( len(self.trajectory), len(self.trajectory) ))
        for i in range(len(self.trajectory)):
            for j in range(i, len(self.trajectory)):
                self.distanceMatrix[i][j] = Utilities.rmsd(self.trajectory[i], self.trajectory[j])
        self.distanceMatrix = self.distanceMatrix + self.distanceMatrix.transpose()
        pyplot.matshow(self.distanceMatrix)
        pyplot.savefig("unsorted_distance_matrix.png", dpi=500)
        pyplot.clf()
        numpy.save("unsorted_distance_matrix.npy", self.distanceMatrix)
        return self.distanceMatrix

    def _kMedoids(self, k = 10, tmax = 100):
        """ clustering function that returns medoids (array) and clusters (dictionary of indices) """
        distanceMatrix = self.distanceMatrix
        # get distanceMatrix dimensions
        m, n = distanceMatrix.shape

        # if more clusters required than structures provided
        if k > n: 
            raise Exception("k parameter is too large.")

        # random initialization of the indices array
        medoids = numpy.arange(n)
        numpy.random.shuffle(medoids)
        medoids = numpy.sort(medoids[:k])

        # new indices array for calculations
        medoids_new = numpy.copy(medoids)

        #dictionary of clusters
        clusters = {}

        for t in xrange(tmax):
            J = numpy.argmin(distanceMatrix[:,medoids], axis = 1)
            for c in range(k):
                clusters[c] = numpy.where(J == c)[0]
            # medoids update
            for c in range(k):
                J = numpy.mean(distanceMatrix[numpy.ix_(clusters[c], clusters[c])], axis = 1)
                j = numpy.argmin(J)
                medoids_new[c] = clusters[c][j]
            numpy.sort(medoids_new)
            # convergence check
            if numpy.array_equal(medoids, medoids_new):
                break
            medoids = numpy.copy(medoids_new)
        else:
            # cluster update
            J = numpy.argmin(distanceMatrix[:,medoids], axis=1)
            for kappa in range(k):
                clusters[kappa] = numpy.where(J==kappa)[0]

        # result: M (medoids) and C (clusters)
        print("Medoids (unsorted): ", medoids)
        print("Clusters (unsorted): ", clusters)
        return medoids, clusters

    def perform_clustering(self):
        self._calculateDistanceMatrix()
        medoids, clusters = self._kMedoids()
        self.medoids_structures = Utilities._reorder(self.trajectory, medoids)
        self.clusters_list = []
        self.labels = []
        for cluster in xrange(max(clusters.keys())):
            cluster_structures = Utilities._reorder(self.trajectory, clusters[cluster])
            self.clusters_list.append(cluster_structures)
            for element in clusters[cluster]:
                self.labels.append(element)
        # reorder the distance matrix
        self.distanceMatrix = self.distanceMatrix[self.labels][:,self.labels]
        numpy.save("sorted_distance_matrix.npy", self.distanceMatrix)
        pyplot.matshow(self.distanceMatrix)
        pyplot.savefig("sorted_distance_matrix.png", dpi=500)
        pyplot.clf()
        #zapisywanie modeli
        number = 1
        for medoid in self.medoids_structures:
            output_file = open("sa_model_"+str(number)+".pdb", 'w')
            output_file.write(medoid.atoms.make_pdb())
            output_file.close()
            number += 1

        #zapisywanie klastrow
        number = 1
        for cluster in self.clusters_list:
            number_2 = 1
            for model in cluster:
                # output_file = open("sa_cluster_"+str(number)+"_"+str(number_2)+".pdb", 'w')
                # output_file.write(model.atoms.make_pdb())
                # output_file.close()
                number_2 += 1
            number += 1

        return self.medoids_structures, self.clusters_list, self.labels

import numpy

from cabsDock.trajectory import Trajectory


class Clustering(object):
    """
    Clustering is a class for performing structural clustering of the models according to similarity of selection.

    If selection is not provided, whole model is used.

    Returns clusters as a list of Trajectory objects with additional 'cluster_quality' as an atribute.

    """

    def __init__(self, trajectory, selection):
        super(Clustering, self).__init__()
        if selection:
            self.trajectory = trajectory.select(selection)
        else:
            self.trajectory = trajectory
        self.distance_matrix = None

    def calculate_distance_matrix(self):
        self.distance_matrix = self.trajectory.rmsd_matrix()
        return self.distance_matrix

    def k_medoids(self, k, tmax=100):
        """
        Performs k-medoid clustering. Retruns medoid indices and a cluster indices  dictionary.

        :param k: number of medoids.
        :param tmax: maximum number of iterations.
        :return: medoid_ndx (list of medoid indices), clusters (dictionary of indeces of cluster elements)
        """
        distance_matrix = self.distance_matrix
        m, n = distance_matrix.shape
        if k > n:
            raise Exception(
                'The number of medoids {0} exceeds the number of structures to be clustered{1}'.format(
                    k, n
                )
            )
        medoid_ndx = numpy.arange(n)
        numpy.random.shuffle(medoid_ndx)
        medoid_ndx = numpy.sort(medoid_ndx[:k])
        medoid_ndx_new = numpy.copy(medoid_ndx)
        clusters = {}
        for t in xrange(tmax):
            j = numpy.argmin(distance_matrix[:, medoid_ndx], axis=1)
            for _k in range(k):
                clusters[_k] = numpy.where(j == _k)[0]
            for _k in range(k):
                j = numpy.mean(distance_matrix[numpy.ix_(clusters[_k], clusters[_k])], axis=1)
                j_min = numpy.argmin(j)
                medoid_ndx_new[_k] = clusters[_k][j_min]
                numpy.sort(medoid_ndx_new)

            if numpy.array_equal(medoid_ndx, medoid_ndx_new):
                break
            medoid_ndx = numpy.copy(medoid_ndx_new)
        else:
            j = numpy.argmin(distance_matrix[:, medoid_ndx], axis=1)
            for _k in range(k):
                clusters[_k] = numpy.where(j == _k)[0]

        # return results
        return medoid_ndx, clusters

    def cabs_clustering(self):
        """
        Performs default cabs-dock clustering (10-medoids) and returns clusters as list of 10 trajectory.Trajectory object with
        trajectory.coordinates.shape = [1, n_models_in_cluster, n_atoms, 3] and medoids as a trajectory.Trajectory object with
        trajectory.coordinates.shape = [1, 10, n_atoms, 3].
        :return:
        """
        self.calculate_distance_matrix()
        medoid_ndx, clusters = self.k_medoids(10)
        model_length = len(self.trajectory.template)
        models = self.trajectory.coordinates.reshape(-1, model_length, 3)
        medoids = Trajectory(
            self.trajectory.template,
            numpy.array([models[medoid_ndx, :, :]]),
            [self.trajectory.headers[i] for i in medoid_ndx]
        )
        clusters_as_clusters = []
        for cluster in clusters.keys():
            this_cluster = Cluster(
                self.trajectory.template,
                numpy.array([models[clusters[cluster], :, :]]),
                [self.trajectory.headers[i] for i in clusters[cluster]]
            )
            clusters_as_clusters.append(this_cluster)
        sorting_ndx =  (sorted(range(len(clusters_as_clusters)), key=lambda x: clusters_as_clusters[x].score, reverse = True))
        medoids.coordinates = medoids.coordinates[:, sorting_ndx, :, :],
        medoids.headers = [medoids.headers[i] for i in sorting_ndx]
        clusters_as_clusters = [clusters_as_clusters[i] for i in sorting_ndx]
        # print([cluster.score for cluster in clusters_as_clusters])
        return medoids, clusters


class Cluster(Trajectory):
    """ docstring for Cluster """

    def __init__(self, template, coordinates, headers):
        super(Cluster, self).__init__(template, coordinates, headers)
        self.score = self.get_score()

    def get_score(self):
        # standard CABSdock cluster scoring based on the cluster density.
        if self.coordinates.shape[1]>1:
            score = self.coordinates.shape[1] / numpy.max(self.rmsd_matrix())
        else:
            score = 0 # One-element clusters are assigned 0 score.
        return score

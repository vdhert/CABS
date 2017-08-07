import numpy

from cabsDock.trajectory import Trajectory


class Clustering(object):
    def __init__(self, trajectory, selection):
        """Clustering is a class for performing structural clustering of the models according to similarity of selected
        atom subset. If selection is not provided, whole model is used.
        :param trajectory: :class:'trajectory.Trajectory' instance
        :param selection: string representing the selection (i.e. 'chain A and chain B')
        """
        super(Clustering, self).__init__()
        self.trajectory = trajectory
        self.selection = selection
        self.distance_matrix = None

    def calculate_distance_matrix(self):
        """
        Calculates the distance matrix. Method for future development of methods based on other distance definitions.
        :return: np.array representing a matrix of distances
        """
        self.distance_matrix = self.trajectory.select(self.selection).rmsd_matrix()
        return self.distance_matrix

    def k_medoids(self, k, tmax=100):
        """
        Performs k-medoid clustering. Returns medoid indices and a cluster indices  dictionary.

        :param k: number of medoids.
        :param tmax: maximum number of iterations.
        :return: list medoid indices, dictionary indeces of cluster elements.
        """
        distance_matrix = self.distance_matrix
        m, n = distance_matrix.shape
        if k > n:
            raise Exception(
                'The number of medoids {0} exceeds the number of structures to be clustered {1}'.format(
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

    def cabs_clustering(self, number_of_medoids, number_of_iterations):
        """
        Performs default cabs-dock clustering (10-medoids) and returns clusters as list of 10 trajectory.
        Trajectory object with trajectory.coordinates.shape = [1, n_models_in_cluster, n_atoms, 3] and medoids as
        a trajectory.Trajectory object with trajectory.coordinates.shape = [1, 10, n_atoms, 3].
        :return: trajectory.Trajectory representing medoids, dictionary representing the cluster indeces, list of
        trajectory.Trajectory instances representing the clusters.
        """
        self.calculate_distance_matrix()
        medoid_ndx, clusters = self.k_medoids(number_of_medoids, tmax=number_of_iterations)
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
        sorting_ndx = (
            sorted(range(len(clusters_as_clusters)), key=lambda x: clusters_as_clusters[x].score, reverse=True))
        # print('sorted clusters')
        # print([clusters_as_clusters[i].score for i in sorting_ndx])
        medoids.coordinates = medoids.coordinates[:, sorting_ndx, :, :]
        medoids.headers = [medoids.headers[i] for i in sorting_ndx]
        clusters_as_clusters = [clusters_as_clusters[i] for i in sorting_ndx]
        return medoids, clusters, clusters_as_clusters


class Cluster(Trajectory):
    def __init__(self, template, coordinates, headers):
        """
        Class for representing a structural cluster.
        :param template:
        :param coordinates:
        :param headers:
        """
        super(Cluster, self).__init__(template, coordinates, headers)
        self.score = self.get_score()

    def get_score(self, method='density'):
        """
        Method for cluster scoring. For future development. Now supports two cluster density definitions, based on
        maximal dissimilarity within the cluster (standard) and standard deviation of the dissimilarity (stdev).
        :param method:
        :return:
        """

        def density(cluster, mode='standard'):
            modes = {
                'standard': numpy.max,
                'stdev': numpy.std,
            }
            return cluster.coordinates.shape[1] / modes[mode](cluster.rmsd_matrix())

        methods = {
            'density': density
        }

        if self.coordinates.shape[1] == 1:
            score = 0
        else:
            score = methods[method](self, mode='standard')
        return score

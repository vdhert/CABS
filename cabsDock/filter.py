import numpy
from trajectory import Trajectory

class Filter(object):
    """
    Class for performing trajectory filtering according to chosen criteria (currently -- lowest energy).

    """
    def __init__(self, trajectory,  N = 1000):
        super(Filter, self).__init__()
        self.trajectory = trajectory
        self.N = N

    def filter(self):
        model_length = len(self.trajectory.template)
        models = self.trajectory.coordinates.reshape(-1, model_length, 3)
        model_energies = [header.get_energy() for header in self.trajectory.headers]
        filtered_ndx = self.mdl_fltr(models, model_energies)
        traj = Trajectory(self.trajectory.template, numpy.array([models[filtered_ndx, :, :]]),
                          [self.trajectory.headers[i] for i in filtered_ndx])
        return traj, filtered_ndx

    @staticmethod
    def mdl_fltr(mdls, enrgs, N=None):
        if N is None:
            N = len(mdls)
        low_energy_ndxs = numpy.argsort(enrgs)
        # print(enrgs)
        # print(numpy.array(enrgs)[low_energy_ndxs])
        if len(mdls) <= N:
            filtered_ndx = low_energy_ndxs
        else:
            filtered_ndx = low_energy_ndxs[:N]
        return filtered_ndx

    def cabs_filter(self):
        n_replicas = self.trajectory.coordinates.shape[0]
        n_models = self.trajectory.coordinates.shape[1]
        fromeach = int(self.N / n_replicas)
        filtered_models = []
        filtered_headers = []
        filtered_total_ndx = []
        
        for i, replica in enumerate(self.trajectory.coordinates):
            energies = [header.get_energy(number_of_peptides=self.trajectory.number_of_peptides) for header in self.trajectory.headers if header.replica == i+1]
            headers = [header for header in self.trajectory.headers if header.replica == i+1]
            filtered_ndx = self.mdl_fltr(replica, energies, N=fromeach)
            if len(filtered_models)==0:
                filtered_models = replica[filtered_ndx, :, :]
                filtered_headers = [headers[k] for k in filtered_ndx]
            else:
                filtered_models = numpy.concatenate([filtered_models, replica[filtered_ndx, :, :]])
                filtered_headers += [headers[k] for k in filtered_ndx]
            filtered_total_ndx.extend(numpy.array(filtered_ndx)+i*n_models)
        traj = Trajectory(self.trajectory.template, numpy.array([filtered_models]), filtered_headers)
        return traj, filtered_total_ndx


# mdls = numpy.zeros(10)
# enrgs = numpy.arange(-10,10)
# print enrgs
# numpy.random.shuffle(enrgs)
# print enrgs
# ndx = Filter.mdl_fltr(mdls,enrgs,N=3)
# print(enrgs[ndx])

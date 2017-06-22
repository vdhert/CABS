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
        low_energy_ndxs = numpy.argsort(model_energies)
        if len(models) <= self.N:
            filtered_ndx = low_energy_ndxs
        else:
            filtered_ndx = low_energy_ndxs[:self.N]
        traj = Trajectory(self.trajectory.template, numpy.array([models[filtered_ndx, :, :]]), [self.trajectory.headers[i] for i in filtered_ndx])
        return traj, filtered_ndx
"""
CABSDock module for contact map analysis of trajectories.

Created on 4 June 2017 by Tymoteusz hert Oleniecki.
"""
import numpy
import operator

import matplotlib.pyplot
import matplotlib.ticker

class ContactMapFactory(object):
    def __init__(self, chain1, chain2, temp):
        """Builder for ContactMap.

        Arguments:
        chain1 -- str; char for 1st chain.
        chain2 -- str; char or 2nd chain.
        temp -- cabsDock.atom.Atoms instance containing both given chains.

        """
        chs = {}
        for n, i in enumerate(temp.atoms):
            chs.setdefault(i.chid, []).append(n)
        self.dims = map(len, [chs[chain1], chs[chain2]])
        self.inds1 = reduce(operator.add, [chs[i] for i in chain1])
        self.inds2 = reduce(operator.add, [chs[i] for i in chain2])
        self.ats1 = [temp.atoms[i] for i in self.inds1]
        self.ats2 = [temp.atoms[i] for i in self.inds2]

    def mk_cmap(self, traj, thr, frames=(0, None), replicas=None):
        """Creates map of contacts between two given chains.

        Arguments:
        traj -- numpy.array of propper shape, i.e. Nreplicas x Nframes x Natoms x 3.
        thr -- float; threshold for side chain distance contact.
        frames -- tuple of ints; range of frames to be taken.
        replicas -- list of replicas' indexes to be taken.

        Returns list of ContactMap for each replica in trajectory.
        """
        if replicas is None:
            replicas = range(traj.shape[0])
        resl = []
        for rep in traj[replicas,]:
            cmtx = self.mk_cmtx(self.mk_dmtx(rep[frames[0]]), thr)
            for fra in rep[frames[0] + 1:frames[1]]:
                ncmtx = self.mk_cmtx(self.mk_dmtx(fra), thr)
                cmtx += ncmtx
            resl.append(ContactMap(cmtx, self.ats1, self.ats2))
        return resl

    def mk_cmtx(self, mtx, thr):
        """Returns boolean numpy.array of contacts from given mtx of distances.

        Arguments:
        mtx -- numpy.array of distances between atoms.
        thr -- thresohld below which atoms are in contact.
        """
        return numpy.clip(numpy.sign(-mtx + thr), 0, 1)

    def mk_dmtx(self, vec):
        """Returns 2D numpy.array of distances between atoms from vector of coords.

        Arguments:
        vec -- slice of trajectory.

        Calculates distances between atoms from given chains.
        """
        mtx = numpy.empty(self.dims, dtype=numpy.float)
        for d1, at1 in enumerate(vec[self.inds1,]):
            for d2, at2 in enumerate(vec[self.inds2,]):
                mtx[d1, d2] = numpy.linalg.norm(at1 - at2)
        return mtx


class ContactMap(object):
    def __init__(self, mtx, atoms1, atoms2):
        """Contact map init.

        Arguments:
        mtx -- 2D numpy.array of distances between (pseudo)atoms.
        temp -- cabsDock.atom.Atoms instance; template for cmap.
        """
        self.cmtx = mtx
        self.s1 = atoms1
        self.s2 = atoms2

    def save(self, fname):
        fig = matplotlib.pyplot.figure(figsize=(20, 20))
        sfig = fig.add_subplot(111)

        sfig.matshow(
            self.cmtx.T,
            cmap=matplotlib.pyplot.cm.Oranges
            )
        sfig.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        sfig.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        sfig.tick_params(axis='both', which='major', labelsize=6)
        for atoms, fx in ((self.s1, sfig.set_xticklabels), (self.s2, sfig.set_yticklabels)):
            fx([''] + [(a.chid + str(a.resnum) + a.icode).strip() for a in atoms])
        matplotlib.pyplot.savefig(fname)

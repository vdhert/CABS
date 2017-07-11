"""
CABSDock module for contact map analysis of trajectories.

Created on 4 June 2017 by Tymoteusz hert Oleniecki.
"""
import numpy
import operator

from cabsDock.utils import _chunk_lst
from cabsDock.utils import _extend_last
from cabsDock.utils import _fmt_res_name
from cabsDock.plots import mk_histos_series

import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.colors


class ContactMapFactory(object):
    def __init__(self, chains1, chains2, temp):
        """Builder for ContactMap.

        Arguments:
        chains1 -- list or str; chars for 1st chain(s).
        chains2 -- list or str; chars for 2nd chain(s).
        temp -- cabsDock.atom.Atoms instance containing both given chains.

        """
        chs = {}
        for n, i in enumerate(temp.atoms):
            chs.setdefault(i.chid, []).append(n)
        self.dims = (sum(map(len, [chs.get(ch1, []) for ch1 in chains1])), sum(map(len, [chs.get(ch2, []) for ch2 in chains2])))
        self.inds1 = reduce(operator.add, [chs.get(i, []) for i in chains1])
        self.inds2 = reduce(operator.add, [chs.get(i, []) for i in chains2])
        self.ats1 = [temp.atoms[i] for i in self.inds1]
        self.ats2 = [temp.atoms[i] for i in self.inds2]

    def mk_cmap(self, traj, thr, frames=None, replicas=None):
        """Creates map of contacts between two given chains.

        Arguments:
        traj -- numpy.array of propper shape, i.e. Nreplicas x Nframes x Natoms x 3.
        thr -- float; threshold for side chain distance contact.
        frames -- tuple of ints; indexes of frames to be taken. All frames are taken by default.
        replicas -- list of replicas' indexes to be taken. All replicas are taken by default.

        Returns list of ContactMap for each replica in trajectory.
        """
        
        if replicas is None:
            replicas = xrange(traj.shape[0])
        if frames is None:
            fstf = 0
            frames = slice(1, None)
        else:
            frames = list(frames)
            fstf = frames.pop(0)
        resl = []
        for rep in traj[replicas,]:
            cmtx = self.mk_cmtx(self.mk_dmtx(rep[fstf]), thr)
            nframes = 1
            for fra in rep[frames,]:
                ncmtx = self.mk_cmtx(self.mk_dmtx(fra), thr)
                cmtx += ncmtx
                nframes += 1
            resl.append(ContactMap(cmtx, map(_fmt_res_name, self.ats1), map(_fmt_res_name, self.ats2), nframes))
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
    def __init__(self, mtx, atoms1, atoms2, n):
        """Contact map init.

        Arguments:
        mtx -- 2D numpy.array of distances between (pseudo)atoms.
        atoms1, atoms2 -- cabsDock.atom.Atoms instance; template for cmap.
        n -- number of frames.
        """
        self.cmtx = mtx
        self.s1 = atoms1
        self.s2 = atoms2
        self.n = n

    def zero_diagonal(self):
        numpy.fill_diagonal(self.cmtx, 0)

    def save_fig(self, fname, fmt='svg', norm_n=False):
        """Saves cmap as matrix plot.

        Arguments:
        fname -- str; file name.
        fmt -- str; by default 'svg'. See matplotlib.pyplot.savefig for more inforamtion.
        norm_n -- bool; if True cmap will be normalized by cmap number of frames.
        """
        #~ grid = matplotlib.pyplot.GridSpec(*(numpy.array([5, 0]) + self.cmtx.shape))
        grid = matplotlib.pyplot.GridSpec(2, 1,
                                    width_ratios=(self.cmtx.shape[1],),
                                    height_ratios=(self.cmtx.shape[0], 1))
        sfig = matplotlib.pyplot.subplot(grid[0, :])
        vmax = self.n if norm_n else numpy.max(self.cmtx)
        colors = matplotlib.colors.LinearSegmentedColormap.from_list('bambi',
                ['#ffffff', '#ff0066', '#33cc33', '#999999', '#ff9933'])

        sfig.matshow(
            self.cmtx.T,
            cmap=colors,
            vmin=0.,
            vmax=vmax,
            )
        sfig.tick_params(axis='both', which='major', labelsize=6)
        for lbls, lab_fx, tck_setter, deg in ((self.s1, sfig.set_xticklabels, sfig.yaxis, 90), (self.s2, sfig.set_yticklabels, sfig.xaxis, 0)):
            mjr_loc = int(round(len(lbls)/50.)*50) / 50
            if mjr_loc == 0: mjr_loc = 1
            tck_setter.set_major_locator(matplotlib.ticker.MultipleLocator(mjr_loc))
            lab_fx([''] + lbls[::mjr_loc], rotation=deg, fontsize=4)

        ax2 = matplotlib.pyplot.subplot(grid[1, 0])
        cb = matplotlib.colorbar.ColorbarBase(ax2,
                                cmap=colors,
                                orientation='horizontal',
                                boundaries=range(int(vmax))
                                )

        matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
        matplotlib.pyplot.close()

    def save_histo(self, fname):
        inds1, inds2 = map(sorted, map(set, numpy.nonzero(self.cmtx)))
        inds1lst = _chunk_lst(inds1, 15)
        trg_vls = [[numpy.sum(self.cmtx[i,:]) for i in inds] for inds in inds1lst]
        vls = [[numpy.sum(self.cmtx[:,i]) for i in inds2]]
        try:
            _extend_last(trg_vls, 15, 0)
        except IndexError:
            trg_vls = [[0.] * 15]
        vls.extend(trg_vls)
        lbls = [[self.s2[i] for i in inds2]]
        trg_lbls = [[self.s1[i] for i in inds] for inds in inds1lst]
        try:
            _extend_last(trg_lbls, 15, "")
        except IndexError:
            trg_lbls = [[""] * 15]
        lbls.extend(trg_lbls)
        mk_histos_series(vls, lbls, fname)

    def save_txt(self, stream):
        """Saves contact list in CSV format.

        Argument:
        stream -- file-like object; stream to which text will be passed.
        """
        inds1, inds2 = numpy.nonzero(self.cmtx)
        stream.write("# n=%i\n" % self.n)
        for m1, m2, (c1, c2) in zip([self.s1[i] for i in inds1], [self.s2[i] for i in inds2], zip(inds1, inds2)):
           stream.write("%s\t%s\t%.3f\n" % (m1, m2, self.cmtx[c1, c2] / float(self.n)))

    def save_all(self, fname):
        """Creates txt and png of given name."""
        with open(fname + '.txt', 'w') as f:
            self.save_txt(f)
        self.save_fig(fname)

    def __add__(self, other):
        """Addition of cmaps sums their matrices.

        Raises ValueError for cmaps of different particles.
        """
        if self.s1 != other.s1 or self.s2 != other.s2:
            raise ValueError("Cannot sum different particles' contact maps.")
        return ContactMap(self.cmtx + other.cmtx, self.s1, self.s2, self.n + other.n)

import matplotlib.pyplot
from matplotlib.ticker import MaxNLocator

import numpy
from itertools import chain

from cabsDock.utils import _chunk_lst


def mk_discrete_plot(splot, xvals, series, xlim=None, ylim=None, joined=False):
    """
    Arguments:
    splot -- matplotlib.pyplot.Axes instance; figure subplot to plot on.
    xvals -- sequence of x-axis values.
    series -- nested sequence of data series (floats or ints) of x-axis values len.
    """
    fmt = "o" + (":" if joined else " ")
    for sser, xsvals in zip(series, xvals):
        splot.plot(xsvals, sser, fmt, markersize=1, linewidth=1)
    if xlim:
        splot.set_xlim(xlim)
    if ylim:
        splot.set_ylim(ylim)
    return splot

def drop_csv_file(fname, columns, fmts="%s"):
    """
    Creates *fname* csv file and writes given columns to this file.

    Arguments:
    fname -- name of file to be created.
    columns -- sequences of data sequences. Lists will be truncated to len of the shortest one.
    fmts -- str or sequence of str. C-style string formats to be used (e.g. "%s", "%.3f", ...).
    If only one string is given -- same format is used for all data.
    Otheriwse subsequent fmts are used for corresponding solumns.
    """
    if type(fmts) is str:
        fmts = [fmts for dummy in columns]
    with open(fname + '.csv', 'w') as f:
        for vals in zip(*columns):
            f.write("\t".join([fmt % val for fmt, val in zip(fmts, vals)]))
            f.write('\n')

def plot_E_RMSD(trajectories, rmsds, fname, fmt='svg'):
    """
    Creates energy(RMSD) plots.

    Arguments:
    trajectories -- sequence of trajectories to be used.
    rmsds -- nested sequence of RMDSs.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See matplotlib.pyplot.savefig for more formats.

    Plots figure with three subplots: total and internal energy vs RMSD
    and histogram of RMSDs. All three will be plotted for given data series,
    so nested arrays or lists are expected. Plots will be written in given format
    to file of given name.
    """
    fig, sfigarr = matplotlib.pyplot.subplots(3)
    for ind, etp in zip((0, 1), ('total','interaction')):
        data = [[h.get_energy(mode=etp, number_of_peptides=traj.number_of_peptides) for h in traj.headers] for traj in trajectories]
        xlim = (0, max(chain((5,), *rmsds)))
        ylim = (min(chain((0,), *data)), max(chain((5,), *data)))
        mk_discrete_plot(sfigarr[ind], rmsds, data, xlim, ylim)
        drop_csv_file(fname + "_%s" % etp, (rmsds[0], data[0]), fmts="%.3f")
    for traj, rmsd_list in zip(trajectories, rmsds):
        sfigarr[2].hist(rmsd_list, int(numpy.max(rmsd_list) - numpy.min(rmsd_list)))
    fig.get_axes()[-1].set_xlabel('RMSD')

    matplotlib.pyplot.savefig(fname + '.'+fmt, format=fmt)
    matplotlib.pyplot.close(fig)

def plot_RMSD_N(rmsds, fname, fmt='svg'):
    """Plots and saves to a file RMSD(Nframe) plot.

    Arguments:
    rmsds -- nested sequence of RMSDs.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See matplotlib.pyplot.savefig for more formats.
    """
    for n, rmsd_lst in enumerate(rmsds):
        tfname = fname + '_replica_%i' % n
        nfs = range(len(rmsd_lst))

        fig, sfig = matplotlib.pyplot.subplots(1)
        mk_discrete_plot(sfig, [nfs], [rmsd_lst], joined=True)
        fig.get_axes()[0].set_ylabel('RMSD')
        fig.get_axes()[0].set_xlabel('frame')
        fig.get_axes()[0].yaxis.set_major_locator(MaxNLocator(integer=True))

        matplotlib.pyplot.savefig(tfname + '.' + fmt, format=fmt)
        matplotlib.pyplot.close(fig)
        drop_csv_file(tfname, (map(str, nfs), rmsd_lst), fmts=("%s", "%.3f"))

def graph_RMSF(trajectory, chains, fname, hist=False):
    if hist:
        rmsf_vals = _chunk_lst(trajectory.rmsf(self.initial_complex.receptor_chains), 15, 0)
        lbls = _chunk_lst([i.chid + str(i.resnum) + i.icode for i in trajectory.template.atoms if i.chid in self.initial_complex.receptor_chains], 15, "")
        mk_histos_series(rmsf_vals, lbls, fname + '_hist')
    else:
        rmsf_vals = [trajectory.rmsf(chains)]
        lbls = [[i.chid + str(i.resnum) + i.icode for i in trajectory.template.atoms if i.chid in chains]]
        plot_RMSF_seq(rmsf_vals, lbls, fname + '_seq')
    drop_csv_file(fname, (tuple(chain(*lbls)), tuple(chain(*rmsf_vals))), fmts=("%s", "%.3f"))

def plot_RMSF_seq(series, labels, fname, fmt='svg'):
    """
    Arguments:
    series -- list of sequences of data to be plotted.
    labels -- corresponding ticks labels.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See matplotlib.pyplot.savefig for more formats.
    """
    fig, sfig = matplotlib.pyplot.subplots(1)
    mk_discrete_plot(sfig, map(range, map(len, series)), series, joined=True)
    sfig.set_ylabel('RMSF')
    sfig.set_xticklabels(labels)
    matplotlib.pyplot.savefig(fname + '.'+fmt, format=fmt)
    matplotlib.pyplot.close(fig)


def mk_histos_series(series, labels, fname, fmt='svg'):
    """
    Arguments:
    series -- list of sequences of data to be plotted.
    labels -- corresponding ticks labels.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See matplotlib.pyplot.savefig for more formats.
    """
    fig, sfigarr = matplotlib.pyplot.subplots(len(series), squeeze=False)

    try:
        ylim = max(chain((5,), *series))
    except ValueError:  # for empty series
        ylim = 5
    get_xloc = lambda x: [.5 * i for i in range(len(x))]

    fig.set_figheight(len(series))

    for n, (vls, ticks) in enumerate(zip(series, labels)):
        xloc = get_xloc(ticks)
        sfigarr[n, 0].bar(xloc, vls, width=.51, color='orange')
        sfigarr[n, 0].set_ylim([0, ylim])
        sfigarr[n, 0].set_yticklabels(["0.00", "%.2f" % ylim], fontsize=6)
        sfigarr[n, 0].set_xticks(xloc)
        sfigarr[n, 0].set_xticklabels(ticks, fontsize=6)

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
    matplotlib.pyplot.close(fig)
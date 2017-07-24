import matplotlib.pyplot
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import IndexFormatter
from matplotlib.ticker import NullFormatter

import numpy
from itertools import chain

from cabsDock.utils import _chunk_lst

matplotlib.pyplot.rcParams['axes.prop_cycle'] = matplotlib.pyplot.cycler(color=['#666666', '#ff4000'])

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
    max_data = 5
    for ind, etp in zip((0, 1), ('total','interaction')):
        grid = matplotlib.pyplot.GridSpec(42, 2)
        plot = matplotlib.pyplot.subplot(grid[:-12, :])
        histo = matplotlib.pyplot.subplot(grid[-10:, :])

        data = [[h.get_energy(mode=etp, number_of_peptides=traj.number_of_peptides) for h in traj.headers] for traj in trajectories]
        xlim = (0, max(chain((max_data,), *rmsds)))
        ylim = (min(chain(*data)), max(chain(*data)))
        mk_discrete_plot(plot, rmsds, data, xlim, ylim)
        drop_csv_file(fname + "_%s" % etp, (rmsds[0], data[0]), fmts="%.3f")

        #TODO histo is the same in both cases, could be done once at the beginning
        for traj, rmsd_list in zip(trajectories, rmsds):
            diff = numpy.max(rmsd_list) - numpy.min(rmsd_list)
            n_bins = numpy.arange(0, max_data, 1) if diff < max_data else numpy.arange(0, int(diff), 1)
            histo.hist(rmsd_list, n_bins)

        plot.xaxis.set_major_locator(MaxNLocator(10))
        plot.xaxis.set_minor_locator(MaxNLocator(20))
        plot.tick_params(labelbottom='off')
        plot.set_ylabel('%s energy' % etp.capitalize())
        histo.set_xlabel('RMSD')
        histo.xaxis.set_major_locator(MaxNLocator(10))
        histo.xaxis.set_minor_locator(MaxNLocator(20))
        histo.set_xlim(xlim)
        matplotlib.pyplot.savefig(fname + '_%s.' % etp + fmt, format=fmt)
        matplotlib.pyplot.close()

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
        sfig.set_ylabel('RMSD')
        sfig.set_xlabel('Frame index')
        sfig.xaxis.set_major_formatter(IndexFormatter(nfs))

        matplotlib.pyplot.savefig(tfname + '.' + fmt, format=fmt)
        matplotlib.pyplot.close(fig)
        drop_csv_file(tfname, (map(str, nfs), rmsd_lst), fmts=("%s", "%.3f"))

def graph_RMSF(trajectory, chains, fname, hist=False):
    if hist:
        rmsf_vals = _chunk_lst(trajectory.rmsf(self.initial_complex.receptor_chains), 15, 0)
        lbls = _chunk_lst([i.fmt() for i in trajectory.template.atoms if i.chid in self.initial_complex.receptor_chains], 15, "")
        mk_histos_series(rmsf_vals, lbls, fname + '_hist')
    else:
        rmsf_vals = [trajectory.rmsf(chains)]
        lbls = [i.fmt() for i in trajectory.template.atoms if i.chid in chains]
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
    xvals = [range(len(series[0]))]
    mk_discrete_plot(sfig, xvals, series, xlim=(min(xvals[0]), max(xvals[0])), joined=True)
    sfig.set_ylabel('RMSF')
    sfig.set_xlabel('Residue index')
    sfig.xaxis.set_major_locator(MaxNLocator(25, integer=True))
    #~ fnc = lambda x, y: labels[int(x)]
    #~ sfig.xaxis.set_major_formatter(FuncFormatter(fnc))
    sfig.xaxis.set_major_formatter(NullFormatter())
    #~ for tick in sfig.get_xticklabels():
        #~ tick.set_rotation(90)
    matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
    matplotlib.pyplot.close(fig)


def mk_histos_series(series, labels, fname, fmt='svg', n_y_ticks=5):
    """
    Arguments:
    series -- list of sequences of data to be plotted.
    labels -- corresponding ticks labels.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    n_y_ticks -- int; maximal number of tickes on y axis.
    See matplotlib.pyplot.savefig for more formats.
    """
    fig, sfigarr = matplotlib.pyplot.subplots(len(series), squeeze=False)

    try:
        ylim = max(chain((n_y_ticks,), *series)) + 1
    except ValueError:  # for empty series
        ylim = n_y_ticks + 1
    get_xloc = lambda x: [.5 * i for i in range(len(x))]

    fig.set_figheight(len(series))

    for n, (vls, ticks) in enumerate(zip(series, labels)):
        xloc = get_xloc(ticks)
        sfigarr[n, 0].set_ylim((0, ylim))
        sfigarr[n, 0].bar(xloc, vls, width=.51)
        sfigarr[n, 0].yaxis.set_major_locator(MaxNLocator(n_y_ticks))
        sfigarr[n, 0].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "%i" % x))
        sfigarr[n, 0].set_xticks(xloc)
        sfigarr[n, 0].set_xticklabels(ticks)
        sfigarr[n, 0].tick_params(labelsize=6)

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
    matplotlib.pyplot.close(fig)
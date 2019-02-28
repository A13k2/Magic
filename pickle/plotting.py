from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
# import seaborn; seaborn.set()
import numpy as np
import pickle
import pandas as pd
# from magic.analysis import plot
import os

def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file
def files_start_with(start, path='.'):
    list_of_files = []
    for file in files(path):
        if file.startswith(start):
            list_of_files.append(file)
    return sorted(list_of_files)
def parse_filename(file):
    num = []
    for lit in file:
        if lit.isnumeric() and not lit=='0':
            num.append(lit)
    return num

def phase_plane_plot(ax, E, I, nu, title=''):
    E = E/1000.
    I = I/1000.
    surf = ax.plot_surface(I, E, nu, cmap=cm.coolwarm)
    ax.set_xlabel('\n\n$I_{ext}$ (kHz)')
    ax.set_ylabel('\n\n$E_{ext}$ (kHz)')
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r'$\hat{\nu}$')
    ax.view_init(elev=14, azim=-127)
    if title:
        ax.set_title(title)
    return ax

def phase_plane_analysis(folder):
    E, I = pickle.load( open(folder+'/E_I.p', 'rb'))
    pickle_files = files_start_with('e_i', path=folder)
    number_of_js = len(pickle_files)
    fig = plt.figure(figsize=(15, 20))
    n = 1
    for file in pickle_files:
        ax1 = fig.add_subplot(number_of_js, 2, n, projection='3d')
        ax2 = fig.add_subplot(number_of_js, 2, n+1, projection='3d')
        n += 2
        title = parse_filename(file)[0]
        e, i = pickle.load( open(folder+'/'+file, 'rb'))
        phase_plane_plot(ax1, E, I, e, title='Excitatory Neurons for $J_{ee}=%d$' % int(title))
        phase_plane_plot(ax2, E, I, i, title='Inhibitory Neurons for $J_{ee}=%d$' % int(title))
    plt.tight_layout()
    return fig

def distance_firing_rate(df, number_of_bins, title):
    duration = df['Time'].values.max() - df['Time'].values.min()
    print(duration)
    grouped = df.groupby([pd.cut(df['Distance from Center'], number_of_bins), pd.cut(df['Time'], number_of_bins)])
    neurons_per_area = df.groupby(pd.cut(df['Distance from Center'], number_of_bins))['Sender'].nunique()
#     Magic to change index from categorical to right side of interval
    tmp = grouped.count().unstack().divide(neurons_per_area, axis=0).multiply(1000./(duration/number_of_bins))['Distance from Center']
    tmp.index = [a.right for a in tmp.index]
    ax = tmp.plot(title=title)
    ax.tight_layout()
    ax.set_ylabel(r'$\hat{\nu}$')
    return ax

def plot_intersections(index, exc_index, e, i, E, I, file = ''):
    max_i = i.max()
    arange = np.arange(0.,max_i,1.)
    n = len(arange)
    fig = plt.figure()
    plt.plot(I[:,exc_index], i[:,exc_index])
    plt.title(r'$\nu_{ex, ext} = ' + str(E[0,exc_index]) + '$ (Hz)')
    plt.xlabel(r'$\nu_{in, ext}$ (Hz)')
    plt.ylabel(r'$\hat{\nu_{in}}$ (Hz)')
    print('r$\\nu_{exc, ext} = '+str(E[0,exc_index])+'$')
    for ind in index:
        plt.plot([I[ind,exc_index]], [i[ind,exc_index]], 'ro')
        print('r$\\nu_{'+str(ind)+', inh, ext} = '+str(I[ind,exc_index])+'$')
    if file:
        plt.savefig(file)

def autocorr(df):
    from scipy.signal import argrelextrema
    bin_center, a = autoCorr(df)
    a = a/a.max()
    # x = bin_center-0.3
    x = bin_center - (bin_center.max()-bin_center.min())/2. - bin_center.min()
    plt.plot(x, a)
    maxs = argrelextrema(a, np.greater)[0]
    mid_max = int(maxs.size/2)
    mins = argrelextrema(a, np.less)[0]
    mid_min = int(mins.size/2)
    plt.plot(x[maxs], a[maxs], 'ro')
    if maxs.size > 2 and mins.size > 2:
        period = x[maxs[mid_max+1]] - x[maxs[mid_max]]
        plt.plot([x[maxs[mid_max]], x[maxs[mid_max+1]]], [a[maxs[mid_max]], a[maxs[mid_max]]], 'g-', label="%.3f" % (period,) )
        periodicity = a[maxs[mid_max+1]] - a[mins[mid_min+1]]
        plt.plot([x[maxs[mid_max+1]], x[maxs[mid_max+1]]], [a[mins[mid_min+1]], a[maxs[mid_max+1]]], 'y-', label="%.3f" % (periodicity))
    plt.xlabel('Time (s)')
    plt.ylabel('ACF')
    plt.legend()

def visualization_times(df, title):
    """
    Visualizes a pandas dataframe for a spike detector. Position is needed.
    Firing Rate for Time interval, plot of 4*3.
    :param df: Pandas dataframe
    :param title: Title of the Plot
    :return: Returns figure and axes from pyplot.subplots()
    """
    df_time_cut = df.groupby(pd.cut(df['Time'], 12))
    fig, axes = plt.subplots(nrows=4, ncols=3, sharex=True, sharey=True)
    fig.suptitle(title)
    # plt.subplots_adjust(top=0.9, hspace=0.25, right=0.84)
    for df_now, ax1 in zip(df_time_cut, axes.flat):
        # ax1.set_title(df_now[0], fontdict={'fontsize': 5})
        ax1.text(0.5, 0.5, df_now[0], va="center", ha="center")
    fig.text(0.5, 0.02, 'X', ha='center')
    fig.text(0.02, 0.5, 'Y', va='center', rotation='vertical')
    # plt.tight_layout()
    return fig, axes

def visualization(df, title):
    """
    Visualizes a pandas dataframe for a spike detector. Position is needed.
    Firing Rate for Time interval, plot of 4*3.
    :param df: Pandas dataframe
    :param title: Title of the Plot
    :return: Returns figure and axes from pyplot.subplots()
    """
    plt.clf()
    df_time_cut = df.groupby(pd.cut(df['Time'], 12))
    fig, axes = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True)
    fig.suptitle(title)
    plt.subplots_adjust(top=0.9, hspace=0.25, right=0.84)
    for df_now, ax1 in zip(df_time_cut, axes.flat):
        x = [i[0] for i in df_now[1]['Position']]
        y = [i[1] for i in df_now[1]['Position']]
        img = ax1.hist2d(x, y, bins=40, cmap=plt.cm.jet, normed=True, vmin=0., vmax=1.)
        ax1.set_title(df_now[0], fontdict={'fontsize': 5})
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(img[3], cax=cbar_ax, )
    fig.text(0.5, 0.02, 'X', ha='center')
    fig.text(0.02, 0.5, 'Y', va='center', rotation='vertical')
    # plt.tight_layout()
    return fig, axes

def visualisation_both(df_ex, df_in, title=None):
    import matplotlib.gridspec as gridspec
    """
    Visualizes two pandas dataframe for a spike detector. Position is needed.
    :param df_ex: Pandas dataframe for excitatory neurons
    :param df_in: Pandas dataframe for inhibitory neurons
    :param title: Title of the Plot
    :return: Return figure
    """
    fig = plt.figure()
    df_ex_time_cut = df_ex.groupby(pd.cut(df_ex['Time'], 12))
    df_in_time_cut = df_in.groupby(pd.cut(df_in['Time'], 12))
    gs0 = gridspec.GridSpec(2, 1)

    gs00 = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
    gs01 = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)

    ax1 = []
    ax2 = []
    for i in range(4):
        for j in range(3):
            ax1.append(plt.subplot(gs00[i,j]))
            ax2.append(plt.subplot(gs01[i,j]))
#     fig.suptitle(title)
    plt.subplots_adjust(top=0.9, hspace=0.25, right=0.84)
    for df_now, ax in zip(df_ex_time_cut, ax1):
        x = [i[0] for i in df_now[1]['Position']]
        y = [i[1] for i in df_now[1]['Position']]
        img = ax.hist2d(x, y, bins=40, cmap=plt.cm.jet, normed=True, vmin=0., vmax=1.)
        ax.set_xticks([])
        ax.set_yticks([])
    for df_now, ax in zip(df_in_time_cut, ax2):
        x = [i[0] for i in df_now[1]['Position']]
        y = [i[1] for i in df_now[1]['Position']]
        img = ax.hist2d(x, y, bins=40, cmap=plt.cm.jet, normed=True, vmin=0., vmax=1.)
        ax.set_xticks([])
        ax.set_yticks([])
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(img[3], cax=cbar_ax, )
    fig.text(0.5, 0.08, 'X', ha='center')
    fig.text(0.48, 0.92, 'Excitatory Neurons', ha='center')
    fig.text(0.48, 0.5, 'Inhibitory Neurons', ha='center')
    fig.text(0.08, 0.5, 'Y', va='center', rotation='vertical')
    return fig

def period(df):
    """
    Calcualtes number of periods for given dataframe by dividing number of zero
    crossings by 2
    :param df: Dataframe to use
    :return float: Number of peaks in df
    """
    y = df[['Sender', 'Time']].values[:,1]
    # create histogram
    bin_heights, bin_borders = np.histogram(y/1000., bins=np.arange(0.1,1.,0.005))
    # calcualte center points from borders of hitogram bins
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    # center around 0
    my_array = bin_heights-bin_heights.mean()
    # calculate zero crossings
    return (((my_array[:-1] * my_array[1:]) < 0).sum())/2.

def variance(df):
    """
    Calculates variance of firing rate from pandas DataFrame that contains
    sender and time of firing neurons
    :param df: Pandas DataFram
    :return float: Variance
    """
    # Get firing times
    y = df[['Sender', 'Time']].values[:,1]
    # Get Rate from binning
    bin_heights, bin_borders = np.histogram(y/1000., bins=np.arange(0.1,1.,0.05))
    # Calcualte number of firing neurons
    number_of_neurons = np.unique(df['Sender'].values).size
    # calculate rate in Hz
    rate_bins = (bin_heights)/(0.05*number_of_neurons)
    return np.var(rate_bins)

def getParams(folder, pickle_file='/parameter.p'):
    """
    Return dict of parameters
    :params folder: Folder wchich contains a pickle file (parameter.p) with the
    parameters inside
    :return dict: parameters
    """
    params = pickle.load( open(folder+pickle_file, 'rb'))
    return params

def printParams(folder):
    """
    Print parameters
    :params folder: Folder wchich contains a pickle file (parameter.p) with the
    parameters inside
    """
    params = getParams(folder)
    for param_key in params.keys():
        print('%s : %s' % (param_key, params[param_key],))

def test_trace_inner(ax, avr_e, avr_i, stim_start=None, stim_end=None, title=None):
    avr_e.plot(ax=ax, label=r'$\hat{\nu_e}$', color='green')
    avr_i.plot(ax=ax, label=r'$\hat{\nu_i}$', color='red')
    ax.set_xlabel(r'Time (ms)')
    ax.set_ylabel(r'$\hat{\nu}$ (Hz)')
    if title:
        ax.set_title(title)
    if stim_start and stim_end:
        if stim_end <= stim_start:
            raise Exception('Stimulation start time is bigger then end time in plotting trace.')
        ax.axvspan(stim_start, stim_end, facecolor='0.2', alpha=0.5)

def trace_plot_all_in_one(files, title='Title', stim_start=None, stim_end=None):
    """
    Plot "Trace" of given DataFrame
    A trace is the avg. firing rate binned in time
    :params files: List with 4 files in following order: [exc inc., exc dec.,
        inh inc., inh dec.]
    :params stim_start: plot stim start and stop
    :params title: Title of the Plot
    """
    def prepare_data(_file):
        avr_e, avr_i = pickle.load(open(_file, 'rb'))
        avr_e.index = [ind.left for ind in avr_e.index]
        avr_i.index = [ind.left for ind in avr_i.index]
        return avr_e, avr_i
    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)
    fig.suptitle(title)
    avr_e, avr_i = prepare_data(files[0])
    test_trace_inner(
            axes[0][0], avr_e, avr_i, stim_start=stim_start, stim_end=stim_end, title='Excitatory increase')
    avr_e, avr_i = prepare_data(files[1])
    test_trace_inner(
            axes[0][1], avr_e, avr_i, stim_start=stim_start, stim_end=stim_end, title='Excitatory decrease')
    avr_e, avr_i = prepare_data(files[2])
    test_trace_inner(
            axes[1][0], avr_e, avr_i, stim_start=stim_start, stim_end=stim_end, title='Inhibitory increase')
    avr_e, avr_i = prepare_data(files[3])
    test_trace_inner(
            axes[1][1], avr_e, avr_i, stim_start=stim_start, stim_end=stim_end, title='Inhibitory decrease')
    plt.tight_layout()
    return fig

def autoCorr(df):
    """
    Returns values for autocorrelation plot of df
    :params df: dataframe
    :return list
    """
    y = df['Time'].values/1000.
    y = y[y > 50./1000.]
    # Get Rate from binning
    binning = np.arange(y.min(),y.max(),0.005)
    bin_heights, bin_borders = np.histogram(y, bins=binning)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    number_of_neurons = np.unique(df['Sender'].values).size
    rate_bins = (bin_heights)/(0.05*number_of_neurons)
    # a = plt.acorr(rate_bins, maxlags=17)
    return bin_centers, np.correlate(rate_bins,rate_bins,'same')

def calculate_periodicity_and_period(df):
    from scipy.signal import argrelextrema
    bin_center, a = autoCorr(df)
    a = a/a.max()
    # x = bin_center-0.3
    x = bin_center - (bin_center.max()-bin_center.min())/2. - bin_center.min()
    maxs = argrelextrema(a, np.greater)[0]
    mins = argrelextrema(a, np.less)[0]
    if maxs.size > 2 and mins.size > 2:
        mid_max = int(maxs.size/2)
        mid_min = int(mins.size/2)
        period = x[maxs[mid_max+1]] - x[maxs[mid_max]]
        periodicity = a[maxs[mid_max+1]] - a[mins[mid_min+1]]
    else:
        period = 0.
        periodicity = 0.
    return periodicity, period

def period_periodicity_plot(folder):
    periodicity_folders = sorted([f[0] for f in os.walk(folder)][1:])
    nus = []
    ex_pys = []
    ex_ps = []
    in_pys = []
    in_ps = []
    for periodicity_folder in periodicity_folders:
        nu = periodicity_folder.split('_')[-1]
        df_ex = pickle.load(open(periodicity_folder+'/df_ex.pickle', 'rb'))
        ex_py, ex_p = calculate_periodicity_and_period(df_ex)
        df_in = pickle.load(open(periodicity_folder+'/df_in.pickle', 'rb'))
        in_py, in_p = calculate_periodicity_and_period(df_in)
        nus.append(nu)
        ex_pys.append(ex_py)
        ex_ps.append(ex_p*1000.)
        in_pys.append(in_py)
        in_ps.append(in_p*1000.)
    nus = [int(nu)/100 for nu in nus]
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    axes[0].plot(nus, ex_pys, label='Excitatory')
    axes[0].plot(nus, in_pys, label='Inhibitory')
    axes[0].set_xlabel(r'$\hat{\nu}_{in, ext}$ (kHz)')
    axes[0].set_ylabel('Periodicity')
    axes[0].legend()
    axes[1].plot(nus[:-3], ex_ps[:-3], label='Excitatory')
    axes[1].plot(nus[:-3], in_ps[:-3], 'r-', label='Inhibitory', dashes=[6,6])
    axes[1].set_xlabel(r'$\hat{\nu}_{in, ext}$ (kHz)')
    axes[1].set_ylabel('Period (ms)')
    axes[1].legend()
    axes[1].set_ylim(bottom=0)
    plt.tight_layout()
    return fig

def get_histo(y, binning, step=0.005, number_of_neurons=16):
    bin_heights, bin_borders = np.histogram(y, bins=binning)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    rate_bins = (bin_heights)/(step*number_of_neurons)
    return bin_centers, rate_bins

def crossCorr(series1, series2):
    min_1 = np.min(series1['Time'].values)
    min_2 = np.min(series2['Time'].values)
    if min_1 <= min_2:
        max_min = min_2
    else:
        max_min = min_1
    series1 = series1[series1['Time'] > max_min]
    series2 = series2[series2['Time'] > max_min]
    y1 = (series1['Time'].values)/1000.
    y2 = (series2['Time'].values)/1000.
    # Get Rate from binning
    binning = np.arange(y1.min(),y1.max(),0.005)
    bin_centers1, rate_bins1 = get_histo(y1, binning)
    bin_centers2, rate_bins2 = get_histo(y2, binning)
    return bin_centers1, np.correlate(rate_bins1,rate_bins2,'same')

def crossCorrPosition(tmp1, tmp2, tmp3):
    fig, ax = plt.subplots()
    bin_center, a = crossCorr(tmp1, tmp2)
    bin_center2, a2 = crossCorr(tmp2, tmp3)
    a = a/a.max()
    a2 = a2/a2.max()
    x_ = bin_center - (bin_center.max()-bin_center.min())/2. - bin_center.min()
    x_2 = bin_center2 - (bin_center2.max()-bin_center2.min())/2. - bin_center2.min()
    ax.plot(x_, a, label="cross(-1,0)")
    ax.plot(x_2, a2, label="cross(0,1)")
    ax.legend()
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('CCF')
    return fig, ax

def crossCorrPopulation(tmp1, tmp2):
    fig, ax = plt.subplots()
    bin_center, a = crossCorr(tmp1, tmp2)
    a = a/a.max()
    x_ = bin_center - (bin_center.max()-bin_center.min())/2. - bin_center.min()
    ax.plot(x_, a)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('CCF')
    return fig, ax

def create_figure_grid(size):
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(20,13))
    gs0 = gridspec.GridSpec(12, 1)
    grid = []
    for i in range(size):
        g0 = gridspec.GridSpecFromSubplotSpec(1, 12, subplot_spec=gs0[2*i], wspace=0.0, hspace=0.0)
        g1 = gridspec.GridSpecFromSubplotSpec(1, 12, subplot_spec=gs0[2*i+1], wspace=0.0, hspace=0.0)
        grid.append([g0,g1])
    return fig, grid

def visualisation_helper(df, g0, print_title=False, nui=3.):
    """
    params: print_title: Boolean, if True print time as title above plots from
                            exc population
    """
    ax1 = []
    ax2 = []
#     for i in range(4):
#         for j in range(3):
#             ax1.append(plt.subplot(g0[0][i,j]))
#             ax2.append(plt.subplot(g0[1][i,j]))
    for j in range(12):
        ax1.append(plt.subplot(g0[0][j]))
        ax2.append(plt.subplot(g0[1][j]))
    df_ex_time_cut = df[0].groupby(pd.cut(df[0]['Time'], 12))
    df_in_time_cut = df[1].groupby(pd.cut(df[1]['Time'], 12))
    q = 0
    for df_now, ax in zip(df_ex_time_cut, ax1):
        x = [i[0] for i in df_now[1]['Position']]
        y = [i[1] for i in df_now[1]['Position']]
        if print_title:
            if q == 0:
                ax.text(-0.77, 0.65, r'Time (ms)', va="center", ha="center")
            ax.set_title(df_now[0], fontdict={'fontsize': 8})
        if q == 0:
            ax.text(-1.3, -0.65, r'$\hat{\nu}_{ext, i}=%.0f$ Hz' % (nui*100., ) , va="center", ha="center")
            ax.text(-0.65, 0., r'Exc.', va="center", ha="center")
        h, x, y, q = ax.hist2d(x, y, bins=40, cmap=plt.cm.jet, normed=True, vmin=0., vmax=1.)
        q.set_edgecolor('face')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    q = 0
    for df_now, ax in zip(df_in_time_cut, ax2):
        if q == 0:
            ax.text(-0.65, 0., r'Inh.', va="center", ha="center")
        x = [i[0] for i in df_now[1]['Position']]
        y = [i[1] for i in df_now[1]['Position']]
        h, x, y, q = ax.hist2d(x, y, bins=40, cmap=plt.cm.jet, normed=True, vmin=0., vmax=1.)
        q.set_edgecolor('face')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    return q

def visualisation_monster(dfs, title=None, nuis=[]):
    import matplotlib.gridspec as gridspec
    fig, grid = create_figure_grid(len(dfs))
    plt.subplots_adjust(top=0.9, hspace=0.25, right=0.84)
    q = 0
    for df, g0, nui in zip(dfs, grid, nuis):
        if q == 0:
            q = visualisation_helper(df, g0, print_title=True, nui = nui)
        else:
            q = visualisation_helper(df, g0, nui = nui)
    return fig

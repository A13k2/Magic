import nest
import nest.topology as tp
import numpy as np
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import matplotlib.animation as animation
import topology_2d_helper as magic


def distanceFiringRate(eventSenders, ctr, min=0, max=0.5, bins=5, neurons_per_gridpoint=8, title='Firing Rate vs. Distance', gridSize=400):
    """
    Plots firing Rates of Time for different Distances
    :param eventSenders: Events perceived from Spike Recording device
    :param ctr: Center Element
    :param min: Minimum distance
    :param max: Maximum distance
    :param bins: Number of bins (distance)
    :param neurons_per_gridpoint: Number of Neurons of the type to be plotted per grid point
    :param title: Title of the Plot
    :param gridSize: Number of grid points in the network
    :return:
    """
    senders = [i for i in eventSenders['senders']]
    times = [i for i in eventSenders['times']]
    distances = [i for i in tp.Distance(senders, [ctr])]
    df = pd.DataFrame({'Sender': senders,
                       'Time': times,
                       'Distance': distances})
    df = df[df['Distance'] <= 0.5]
    grouped = df.groupby([pd.cut(df['Distance'], 10), pd.cut(df['Time'], 10)])
    grouped = grouped.count()['Sender'].unstack()
    #Calculate Neurons per area
    dist = pd.DataFrame({'Distance': [i for i in tp.Distance([i for i in range(2, gridSize+2)], [ctr])]})
    dist = dist[dist['Distance'] <= 0.5]
    neuronsPerArea = dist.groupby(pd.cut(dist['Distance'], 10)).count()['Distance'].multiply(neurons_per_gridpoint)
    grouped = grouped.divide(neuronsPerArea, axis=0)
    grouped.plot(title=title)


def makePandas(senderEvents, ctr):
    """
    Return a pandas DataFrame (Sender, Time, Position, Distance from ctr) for the recordings of a spike detector
    :param senderEvents: Recordings from spike detector
    :param ctr: Center Element(GID)
    :return: Pandas dataframe with Sender, Time, Position and Distance
    """
    senders = [int(i) for i in senderEvents['senders']]  #int necessary, because GID's for tp.GetPosition. Otherwise error
    times = [i for i in senderEvents['times']]
    distances = [i for i in tp.Distance(senders, [ctr])]
    positions = [i for i in tp.GetPosition(senders)]
    return pd.DataFrame({'Sender': senders,
                         'Time': times, 'Distance from Center': distances,
                         'Position': positions
                         })


def convertTopologyToRealID(pandasFrame, gridSize=400):
    for index, row in pandasFrame.iterrows():
        row['Sender'] = row['Sender'] % gridSize
    return pandasFrame


def interspikeIntervals(spikeTrain):
    """
    Returns Interspike Intervals (Time difference between Spikes)
    :param spikeTrain:
    :return:
    """
    return np.diff(spikeTrain)


def fanoFactor(df, bins=10, gridSize=400, tMin=0., tMax=250., tStep=25.):
    """
    Creates Fano Factor Plots (for different times)
    :param df: Pandas dataframe with Sender and Time
    :param bins: Binning in time
    :param gridSize: Number of grid points in the network
    :param tMin: time to start from
    :param tMax: time to end at
    :param tStep: Step
    :return: ls: Fano Factor
    :return: ls: Time for Fano Factor
    """
    df['Sender'] = df['Sender'] % gridSize
    ls = []
    tReturn = []
    t_old = 0
    for t in np.arange(tMin, tMax, tStep):
        s = []
        for i in range(2, gridSize+2):
            oneSender = df[df['Sender'] == i]['Time']
            oneSender = oneSender[oneSender >= t_old]
            oneSender = oneSender[oneSender < t]
            interSpike = np.diff(oneSender)
            if (len(interSpike) > 1):
                fano = np.var(interSpike)/np.mean(interSpike)
                s.append(fano)
        if (s != []):
            s = np.asarray(s)
            ls.append(np.mean(s)/100.)
            tReturn.append(t)
        t_old = t
    return ls, tReturn


def spike_count_histogram(df, t_start, t_stop, t_step, number_of_neurons, grid_size):
    """
    Returns Histogram splitted for different Times
    :param df:
    :param t_start:
    :param t_stop:
    :param t_step:
    :param number_of_neurons:
    :param grid_size:
    :return:
    """
    y, x = np.histogram(df['Time']/1000., np.arange(t_start, t_stop, t_step))
    y = y/((t_stop-t_start)*grid_size*number_of_neurons)
    return x, y


def spike_count_histogram_plot(df, t_start, t_stop, t_step, number_of_neurons, grid_size):
    """
    Return Plot for spike_count_histogram
    :param df:
    :param t_start:
    :param t_stop:
    :param t_step:
    :param number_of_neurons:
    :param grid_size:
    :return:
    """
    x,y = spike_count_histogram(df, t_start, t_stop, t_step, number_of_neurons, grid_size)
    return plt.bar(x[:-1], y, width=0.9*t_step)


def fanoFactorNew(df, t_start, t_stop, t_step, number_of_neurons, grid_size):
    """
    Return Fano Factor for pandas dataframe and for different times
    :param df:
    :param t_start:
    :param t_stop:
    :param t_step:
    :param number_of_neurons:
    :param grid_size:
    :return:
    """
    x, y = spike_count_histogram(df, t_start, t_stop, t_step, number_of_neurons, grid_size)
    return np.var(y)/np.mean(y)

def fanoFactorTime(df, t_start, t_stop, t_step, number_of_neurons, grid_size, bins=10):
    """
    Fano Factor at different Time
    :param df:
    :param t_start:
    :param t_stop:
    :param t_step:
    :param number_of_neurons:
    :param grid_size:
    :param bins:
    :return:
    """
    df_time_cut = df.groupby(pd.cut(df['Time'], bins))
    ts = []
    fano = []
    for t, df_curr in df_time_cut:
        fano_curr = fanoFactorNew(df_curr, t._repr_base()[0]/1000., t._repr_base()[1]/1000., t_step, number_of_neurons, grid_size)
        fano.append(fano_curr)
        ts.append(t._repr_base()[0])
    return ts, fano


def fanoFactorTimePlot(df, t_start, t_stop, t_step, number_of_neurons, grid_size, bins=10):
    """
    Return Plot for fanoFactorTime
    :param df:
    :param t_start:
    :param t_stop:
    :param t_step:
    :param number_of_neurons:
    :param grid_size:
    :param bins:
    :return:
    """
    t_fano, fano = fanoFactorTime(df, t_start, t_stop, t_step, number_of_neurons, grid_size, bins=bins)
    fig = plt.plot(t_fano, fano)
    plt.ylim(0., 1.)
    return fig


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
    return fig, axes


def raster_plot(senders=0, timeS=0, eventSenders=None, timeStamp=0, neurons=0, title=None, gridSize=400):
    """
    Returns Rasterplot of given events.

    ************************
    Parameters:
        eventSenders: Neuron ID & Time of spike
            Get by : nest.GetStatus(....)[0]
        timeStamp: Range of time to be plotted
        neurons: ID's of Neurons to be plotted
        title: Title of the Figure
    ************************
    Returns:
        Nest Raster Plot
    """
    if eventSenders != None:
        senders = eventSenders['senders']
        timeS = eventSenders['times']
    event = []
    ts = []
    for events, times in zip(senders, timeS):
        event.append(events % gridSize)
        ts.append(times)
    if timeStamp == 0:
        timeStamp = ts.copy()
    if neurons == 0:
        neurons = event.copy()
    return nest.raster_plot._make_plot(ts, timeStamp, event, neurons, title=title)


def distance(eventSenders, distanceMin, distanceMax, distanceFrom):
    """
    Returns ID of Sender (Firing Neuron), and related firing time for Neuron in given distance from 'distanceFrom'
    :param eventSenders: Events perceived from Spike Recording device
    :param distanceMin: Minimum distance
    :param distanceMax: Maximum distance
    :param distanceFrom: Defines from wich elemnt to take the distance from
    :return:
        event:  Event ID's
        ts:     Times of spiking
    """
    event = []
    ts = []
    senders = eventSenders['senders']
    times = eventSenders['times']
    for i in np.arange(0, len(senders), 1):
        currSender = senders[i]
        currTime = times[i]
        distance = tp.Distance([currSender], [distanceFrom])[0]
        if (distance >= distanceMin and distance <= distanceMax):
            event.append(currSender)
            ts.append(currTime)
    return event, ts


def recordElectrode(df, posX, posY, numberOfNeurons=1):
    """
    Shows firing Rate for Neuron at given Position
    :return: Plot
    """
    df = df[df['Position'] == (posX, posY)]
    df_grouped_sender = df.groupby('Sender')
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.suptitle('Recorded Spikes at Position: $(' + str(round(posX, 2)) + ', ' + str(round(posY, 2)) + ')$')
    fig.text(0.5, 0.04, '$s$', ha='center')
    fig.text(0.04, 0.5, 'Spikes', va='center', rotation='vertical')
    for (sender, df_now), ax in zip(df_grouped_sender, axes.flat):
        df_curr = pd.DataFrame(df_now)
        df_curr = df_curr.groupby(pd.cut(df_curr['Time'], 40)).count()
        t = [round(t._repr_base()[1]/1000., 4) for t in df_curr.index]
        ax.bar(t, df_curr['Time'], align='center', width=0.05)
#        ax.set_xlabel('$s$')
#        ax.set_ylabel('Spikes')
    return fig, axes


class RandomBalancedNetwork:
    def __init__(self, parameters):
        self.parameters = parameters
        self.gridSize = parameters['Columns']*parameters['Rows']
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": 0.1, "print_time": True, "overwrite_files": True})
        nest.SetKernelStatus({"local_num_threads": 4})
        nest.CopyModel('iaf_psc_alpha', 'exci')
        nest.CopyModel('iaf_psc_alpha', 'inhi')
        nest.CopyModel('static_synapse', 'exc', {'weight': self.parameters['Excitational Weight']})
        nest.CopyModel('static_synapse', 'inh', {'weight': self.parameters['Inhibitory Weight']})
        nest.CopyModel('static_synapse', 'inh_strong', {'weight': self.parameters['Weight Stimulus']})
        self.l = tp.CreateLayer({'rows': self.parameters['Rows'],
                                 'columns': self.parameters['Columns'],
                                 'elements': ['exci', self.parameters['Number excitational cells'], 'inhi', self.parameters['Number inhibitory cells']]})
        cdict_e2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius excitational']}},
                     'kernel': {'gaussian': {'p_center': 0.8, 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': 2.0, 'a': 0.02}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'exc'}
        cdict_e2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Excitational Weight']}},
                     'kernel': {'gaussian': {'p_center': 0.8, 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': 2.0, 'a': 0.02}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'exc'}
        cdict_i2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': 0.8, 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': 2.0, 'a': 0.02}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'inh'}
        cdict_i2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': 0.8, 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': 2.0, 'a': 0.02}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'inh'}
        tp.ConnectLayers(self.l, self.l, cdict_e2i)
        tp.ConnectLayers(self.l, self.l, cdict_e2e)
        tp.ConnectLayers(self.l, self.l, cdict_i2i)
        tp.ConnectLayers(self.l, self.l, cdict_i2e)
        stim = tp.CreateLayer({'rows': 1,
                               'columns': 1,
                               'elements': 'poisson_generator'})
        stim_i = nest.GetLeaves(stim, local_only=True)[0]
        stim_i = nest.GetLeaves(stim, local_only=True)[0]
        nest.SetStatus(stim_i, {'rate': parameters['Background rate']})
        background_stim_dict = {'connection_type': 'divergent',
                                'mask': {'circular': {'radius': 2.}},
                                'synapse_model': 'exc'}
        tp.ConnectLayers(stim, self.l, background_stim_dict)
        stim2 = tp.CreateLayer({'rows': 1,
                                'columns': 1,
                                'elements': 'poisson_generator'})
        self.stim2_i = nest.GetLeaves(stim2, local_only=True)[0]
        nest.SetStatus(self.stim2_i, {'rate': 0.0})
        self.cdict_stim2 = {'connection_type': 'divergent',
                            'kernel': {'gaussian': {'p_center': 1., 'sigma': self.parameters['Sigma Stimulus']}},
                            'mask': {'circular': {'radius': self.parameters['Radius stimulus']},
                                     'anchor': [0., 0.]},
                            'targets': {'model': 'exci'},
                            'synapse_model': 'inh_strong'}
        tp.ConnectLayers(stim2, self.l, self.cdict_stim2)
        rec = nest.Create("spike_detector")
        nrns = nest.GetLeaves(self.l, local_only=True)[0]
        nest.Connect(nrns, rec)
        self.rec_ex = tp.CreateLayer({'rows': 1,
                                      'columns': 1,
                                      'elements': 'spike_detector'})
        cdict_rec_ex = {'connection_type': 'convergent',
                        'sources': {'model': "exci"}}
        tp.ConnectLayers(self.l, self.rec_ex, cdict_rec_ex)
        self.rec_in = tp.CreateLayer({'rows': 1,
                                      'columns': 1,
                                      'elements': 'spike_detector'})
        cdict_rec_in = {'connection_type': 'convergent',
                        'sources': {'model': 'inhi'}}
        tp.ConnectLayers(self.l, self.rec_in, cdict_rec_in)

    def start_simulation(self):
        nest.Simulate(self.parameters['Time before stimulation'])
        nest.SetStatus(self.stim2_i, {'rate': self.parameters['Stimulus rate']})
        nest.Simulate(self.parameters['Time of stimulation'])
        nest.SetStatus(self.stim2_i, {'rate': 0.0})
        nest.Simulate(self.parameters['Time after Stimulation'])
        rec_ex_true = nest.GetLeaves(self.rec_ex, local_only=True)[0]
        rec_in_true = nest.GetLeaves(self.rec_in, local_only=True)[0]
        self.events_ex = nest.GetStatus(rec_ex_true, "events")[0]
        self.events_in = nest.GetStatus(rec_in_true, "events")[0]
        self.df_ex = magic.makePandas(self.events_ex, tp.FindCenterElement(self.l)[0])
        self.df_in = magic.makePandas(self.events_in, tp.FindCenterElement(self.l)[0])

    def writeParametersToFile(self, file):
        """
        Writing Parameters to file
        __________________________
        Writing crucial Parameters from Simulation to an File.
        This makes it easier to understand where the Plots came from.
        __________________________
        Parameters
        parameter: Dictionary of Parameters
        file:      Where to save the file
        """
        with open(file, 'w') as f:
            for para in self.parameters:
                f.write(para+'\t'+str(self.parameters[para])+'\n')
            f.close()
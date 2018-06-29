import nest
import nest.topology as tp
import numpy as np
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import matplotlib.animation as animation
import topology_2d_helper as magic


def writeParametersToFile(parameter, file):
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
    f = open(file, 'w')
    for para in parameter:
        f.write(para+'\t'+str(parameter[para])+'\n')
    f.close()


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
                         'Time': times,
                         'Distance from Center': distances,
                         'Position': positions
                         })

def convertTopologyToRealID(pandasFrame, gridSize=400):
    for index, row in pandasFrame.iterrows():
        row['Sender'] = row['Sender'] % gridSize
    return pandasFrame

def interspikeIntervals(spikeTrain):
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
        img = ax1.hist2d(x, y, bins=40, normed=True)
        ax1.set_title(df_now[0], fontdict={'fontsize': 5})
    #         fig.colorbar(img, ax=ax1)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(img[3], cax=cbar_ax, )
    return fig, axes

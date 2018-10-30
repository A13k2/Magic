import nest
import os
import nest.topology as tp
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import topology_2d_helper as magic

def linear(x,a,b):
    return a+x*b

def gaussian(x, mu, sig, maximum=1.):
    return maximum*np.exp(-np.power(x - mu, 2.) / (2. * np.power(sig, 2.)))


def distancePlotsLazy(Start, End, Step, network, folder):
    for events, title, neurons_folder in zip([network.events_ex, network.events_in], ["Excitatory", "Inhibitory"], ['/excitatory_neurons', '/inhibitory_neurons']):
        for start, end in zip(np.arange(Start, End, Step), np.arange(Start + Step, End + Step, Step)):
            plt.close()
            currentTitle = title + 'Neurons from ' + str(start) + ' to ' + str(end)
            print(currentTitle)
            times, neurons = magic.distance(events, start, end, tp.FindCenterElement(network.l)[0])
            magic.raster_plot(senders=times, timeS=neurons, title=currentTitle, gridSize=network.gridSize)
            plt.savefig(folder + neurons_folder + '/dist_' + str(end) + '_raster_' + title + '.png',
                dpi=300)


def stimulationControlLazy(network, folder):
    for df, title, neurons_folder in zip([network.df_ex, network.df_in], ['Excitatory', 'Inhibitory'], ['/excitatory_neurons', '/inhibitory_neurons']):
        df_before_stim = df[df['Time'] <= network.parameters['Time before stimulation']]
        plt.close()
        magic.visualization(df_before_stim, title + ' Neurons before inhibition, ' + network.parameters['Name'])
        plt.savefig(folder + neurons_folder + '/visu_before_stim_' + title + '_neurons.png',
                    dpi=300)
        df_while_stim = df[df['Time'] > network.parameters['Time before stimulation']]
        df_while_stim = df_while_stim[
            df_while_stim['Time'] <= network.parameters['Time before stimulation'] + network.parameters[
                'Time of stimulation']]
        plt.close()
        magic.visualization(df_while_stim, title + ' Neurons while inhibition, ' + network.parameters['Name'])
        plt.savefig(folder + neurons_folder + '/visu_while_stim_' + title + '_neurons.png',
                    dpi=300)
        df_after_stim = df[
            df['Time'] > network.parameters['Time before stimulation'] + network.parameters['Time of stimulation']]
        plt.close()
        magic.visualization(df_after_stim, title + ' Neurons after inhibition, ' + network.parameters['Name'])
        plt.savefig(folder + neurons_folder + '/visu_after_stim_' + title + '_neurons.png',
                    dpi=300)


def fanoFactorTimeLazy(network, folder):
    plt.clf()
    magic.fanoFactorTimePlot(network.df_ex, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number excitational cells'], network.gridSize, bins=20)
    plt.savefig(folder + '/excitatory_neurons/fano_factor_time.png',
                dpi=300)
    plt.close()
    plt.clf()
    magic.fanoFactorTimePlot(network.df_in, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number inhibitory cells'], network.gridSize, bins=20)
    plt.savefig(folder + '/inhibitory_neurons/fano_factor_time.png',
                dpi=300)
    plt.close()


def spikeCountHistogramLazy(network, folder):
    plt.close()
    magic.spike_count_histogram_plot(network.df_ex, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number excitational cells'], network.gridSize)
    plt.savefig(folder + '/excitatory_neurons/spike_count_histogram.png',
                dpi=300)
    plt.close()
    magic.spike_count_histogram_plot(network.df_in, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number inhibitory cells'], network.gridSize)
    plt.savefig(folder + '/inhibitory_neurons/spike_count_histogram.png',
                dpi=300)


def recordElectrodeLazy(network, folder):
    neurons_x_position_distance = 0.5/(network.parameters['Columns']/2.)
    neurons_x_position = [(i, 0.) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)]
#Add y-Positions
    neurons_x_position += ([(0., i) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)])
#Add diagonal
    neurons_x_position += ([(i, i) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)])
    neurons_x_position = tp.FindNearestElement(network.l, neurons_x_position)
    neurons_position = tp.GetPosition(neurons_x_position)
    for position in neurons_position:
        plt.close()
        magic.recordElectrode(network.df_ex, position[0], position[1], 8)
        plt.savefig(folder + '/excitatory_neurons/electrode_' + str(round(position[0], 4)) + ', ' + str(round(position[1], 4)) + '.png',
                    dpi=300)
        plt.close()
        magic.recordElectrode(network.df_in, position[0], position[1], 8)
        plt.savefig(folder + '/inhibitory_neurons/electrode_' + str(round(position[0], 4)) + ', ' + str(round(position[1], 4)) + '.png',
                    dpi=300)


def recordElectrodeEnviromentLazy(network, folder):
    neurons_x_position_distance = 0.5/(network.parameters['Columns']/2.)
    neurons_x_position = [(i, 0.) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)]
    #Add y-Positions
    neurons_x_position += ([(0., i) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)])
    #Add diagonal
    neurons_x_position += ([(i, i) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)])
    neurons_x_position = tp.FindNearestElement(network.l, neurons_x_position)
    neurons_position = tp.GetPosition(neurons_x_position)
    for position in neurons_position:
        plt.close()
        magic.recordElectrodeEnviroment(network.df_ex, position[0], position[1], 0.05, 0.05)
        plt.savefig(folder + '/excitatory_neurons/electrode_enviroment' + str(round(position[0], 4)) + ', ' + str(round(position[1], 4)) + '.png',
                    dpi=300)
        plt.close()
        magic.recordElectrodeEnviroment(network.df_in, position[0], position[1], 0.05, 0.05)
        plt.savefig(folder + '/inhibitory_neurons/electrode_enviroment' + str(round(position[0], 4)) + ', ' + str(round(position[1], 4)) + '.png',
                    dpi=300)


def rasterPlotLazy(network, folder):
    """
    Creates Raster Plots for excitatory an inhibitory neurons of the given network.
    :param network:
    :param folder:
    :return:
    """
    magic.raster_plot(eventSenders=network.events_ex, gridSize=network.gridSize)
    plt.savefig(folder + '/excitatory_neurons/raster.png', dpi=300)
    magic.raster_plot(eventSenders=network.events_in, gridSize=network.gridSize)
    plt.savefig(folder + '/inhibitory_neurons/raster.png', dpi=300)


def makeDir(folder):
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


def simulationAndAnalysis(parameters, curr_folder='.'):
                delay_visualisation_linear(parameters, curr_folder+'/delay.png')
                weightVisualisation(parameters, curr_folder+'/weights.pdf')
                simulation = magic.RandomBalancedNetwork(parameters)
                simulation.start_simulation()
                simulation.writeParametersToFile(curr_folder + '/parameters.txt')
                makeDir(curr_folder+'/excitatory_neurons')
                makeDir(curr_folder+'/inhibitory_neurons')
                stimulationControlLazy(simulation, curr_folder)
                clusteringPlotLazy(simulation, parameters, curr_folder)
                fanoFactorTimeLazy(simulation, curr_folder)
                recordElectrodeEnviromentLazy(simulation, curr_folder)
                spikeCountHistogramLazy(simulation, curr_folder)
                rasterPlotLazy(simulation, curr_folder)
                # distancePlotsLazy(1000., 2000., 500., simulation, curr_folder)


def weightVisualisation(parameters, file):
    plt.clf()
    x = np.arange(-0.5,0.51,0.01)
    plt.plot(x, gaussian(x, 0., parameters['Sigma excitational'], maximum=parameters['Excitatory Weight Maximum']), label='Excitational')
    plt.plot(x, gaussian(x, 0., parameters['Sigma inhibitory'], maximum=parameters['Inhibitory Weight Maximum']), label='Inhibitory')
    plt.plot(x, gaussian(x, 0., parameters['Sigma excitational'], maximum=parameters['Excitatory Weight Maximum'])-gaussian(x, 0., parameters['Sigma inhibitory'], maximum=parameters['Inhibitory Weight Maximum']), label='Excitational - Inhibitory')
    plt.title('Weight Distribution')
    plt.legend()
    plt.xlabel('Distance from Neuron')
    plt.ylabel('Probability of Connection')
    plt.savefig(file)

def delay_visualisation_linear(parameters, file):
    plt.clf()
    x = np.arange(-0.,0.51,0.01)
    plt.figure(1)
    ax1 = plt.subplot(211)
    plt.title('Delay growth')
    plt.plot(x, linear(x, parameters['e2e delay'],
                       parameters['delay growth multiplier']*parameters['e2e delay']), label='e2e')
    plt.plot(x, linear(x, parameters['e2i delay'],
                       parameters['delay growth multiplier']*parameters['e2i delay']), label='e2i')
    plt.ylabel('Delay in ms')
    plt.legend()
    ax2 = plt.subplot(212, sharex=ax1)
    plt.plot(x, linear(x, parameters['i2i delay'],
                       parameters['delay growth multiplier']*parameters['i2i delay']), label='i2i')
    plt.plot(x, linear(x, parameters['i2e delay'],
                       parameters['delay growth multiplier']*parameters['i2e delay']), label='i2e')
    plt.ylabel('Delay in ms')
    plt.xlabel('Distance from Neuron')
    plt.legend()
    plt.savefig(file)



def clusteringPlotLazy(network, parameters, folder):
    """
    Plots clustering
    """
    def worker1(percentage, folder):
        increased, average, decreased = clusteringPlot(network.df_ex, parameters['Time before stimulation'], 1000., 1500., percChange=percentage)
        exc_folder_during = folder+'/excitatory_neurons/barplots/during/'
        makeDir(exc_folder_during)
        exc_folder_after = folder+'/excitatory_neurons/barplots/after/'
        makeDir(exc_folder_after)
        combinedBarPlot(average, increased, decreased, exc_folder_during+str(percentage)+'.png')
        increased, average, decreased = clusteringPlot(network.df_in, parameters['Time before stimulation'], 1000., 1500., percChange=percentage)
        inh_folder_during = folder+'/inhibitory_neurons/barPlots/during/'
        makeDir(inh_folder_during)
        inh_folder_after = folder+'/inhibitory_neurons/barPlots/after/'
        makeDir(inh_folder_after)
        combinedBarPlot(average, increased, decreased, inh_folder_during+str(percentage)+'.png')
        increased, average, decreased = clusteringPlot(network.df_ex, parameters['Time before stimulation'], 1500., 1700., percChange=percentage)
        combinedBarPlot(average, increased, decreased, exc_folder_after+str(percentage)+'.png')
        increased, average, decreased = clusteringPlot(network.df_in, parameters['Time before stimulation'], 1500., 1700., percChange=percentage)
        combinedBarPlot(average, increased, decreased, inh_folder_after+str(percentage)+'.png')
        return
    for percentage in [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]:
        # t = threading.Thread(target=worker1, args=(percentage, file))
        # t.start()
        worker1(percentage, folder)


def combinedBarPlot(average, increased, decreased, file):
    plt.clf()
    x = np.arange(3)
    all0 = average[0] + increased[0] + decreased[0]
    all1 = average[1] + increased[1] + decreased[1]
    all2 = average[2] + increased[2] + decreased[2]
    data0 = np.array([average[0]/all0*100., average[1]/all1*100., average[2]/all2*100.])
    data1 = np.array([increased[0]/all0*100., increased[1]/all1*100., increased[2]/all2*100.])
    data2 = np.array([decreased[0]/all0*100., decreased[1]/all1*100., decreased[2]/all2*100.])
    p0 = plt.bar(x, data0)
    p1 = plt.bar(x, data1, bottom=data0)
    p2 = plt.bar(x, data2, bottom=data0+data1)
    plt.xticks(x, ('loc.', 'interm.', 'dist.'))
    plt.ylabel('% of all')
    plt.legend((p0[0], p1[0], p2[0]), ('Average', 'Increased', 'Decreased'))
    plt.savefig(file)



def barplotCluster(average, increased, decreased, file="", save=False):
    plt.clf()
    x = np.arange(3)
    all = average + increased + decreased
    data = np.array([average, increased, decreased])
    data = (data/all)*100.
    p = plt.bar(x, data)
    plt.xticks(x, ('avrg.', 'incr.', 'decr.'))
    plt.ylabel(r'% of all ( ' + str(all) + r' )')
    if save and file != "":
        plt.savefig(file)
    return plt


def clusteringPlot(df, warmUpTime, startRecord, stopRecord, percChange=0.05, startwarmUpTime=0.):
    """
    Makes "cluster Plots"
    :param df:
    :param warmUpTime: Time until stimulus/inhibitin is started
    :param startRecord: Time to start recording (doesnt include warmup time)
    :param stopRecord: Time to stop recording
    :param percChange: percentage of change compared to average to classify increased & decreased
    :param startwarmUpTime: Not yet implemented
    :return:
    """
    def calculateAverageFiringRate(df):
        numberNeurons = len(np.unique(df['Sender']))
        return len(df)/float(numberNeurons*(warmUpTime/1000.))
    def splitDataFrame(df):
        df_local = df[df['Distance from Center'] < 0.1]
        df_interm = df[df['Distance from Center'] > 0.1]
        df_interm = df_interm[df_interm['Distance from Center'] < 0.25]
        df_distal = df[df['Distance from Center'] > 0.25]
        df_distal = df_distal[df_distal['Distance from Center'] < 0.45]
        return [df_local, df_interm, df_distal]
    df_warmUp = df[df['Time'] < warmUpTime]
    averageFiringRate = calculateAverageFiringRate(df_warmUp)
    df_interest = df[df['Time'] > startRecord]
    df_interest = df_interest[df_interest['Time'] < stopRecord]
    avr_array = []
    inc_array = []
    dec_array = []
    for df_current in splitDataFrame(df_interest):
    # for df_current in splitDataFrame(df_afterStim):
        average = 0
        decreased = 0
        increased = 0
        for neuron_id in np.unique(df_current['Sender']):
            firingRate = len(df_current[df_current['Sender'] == neuron_id])/((stopRecord-startRecord)/1000.)
            # firingRate = len(df_current[df_current['Sender'] == neuron_id])/(time/1000.)
            firingDifference = firingRate/averageFiringRate
            if (firingDifference > 1.+percChange):
                increased += 1
            elif(firingDifference < 1.-percChange):
                decreased += 1
            else:
                average += 1
        print("Average firing rate:", averageFiringRate)
        all_curr = float(average+decreased+increased)
        print("Number of Neurons with normal Firing Rate: ", average, " (", float(average)/all_curr, ", ", average , "/", all_curr, ")")
        print("Number of Neurons with decreased Firing Rate: ", decreased, " (", float(decreased)/all_curr, ", ", decreased , "/", all_curr, ")")
        print("Number of Neurons with increased Firing Rate: ", increased, " (", float(increased)/all_curr, ", ", increased , "/", all_curr, ")")
        avr_array.append(average)
        inc_array.append(increased)
        dec_array.append(decreased)
    return inc_array, avr_array, dec_array




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
#        import pdb; pdb.set_trace()
        fano_curr = fanoFactorNew(df_curr, t.left/1000., t.right/1000., t_step, number_of_neurons, grid_size)
        fano.append(fano_curr)
        ts.append(t.left)
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


def recordElectrode(*df, posX, posY, numberOfNeurons=1):
    """
    Shows firing Rate for Neuron at given Position
    :return: Plot
    """
    df_pos = df[df['Position'] == (posX, posY)]
    df_grouped_sender = df_pos.groupby('Sender')
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    fig.suptitle('Recorded Spikes at Position: $(' + str(round(posX, 2)) + ', ' + str(round(posY, 2)) + ')$')
    fig.text(0.5, 0.04, '$s$', ha='center')
    fig.text(0.04, 0.5, 'Spikes', va='center', rotation='vertical')
    for (sender, df_now), ax in zip(df_grouped_sender, axes.flat):
        df_curr = pd.DataFrame(df_now)
        df_curr = df_curr.groupby(pd.cut(df_curr['Time'], 40)).count()
        t = [round(t.right/1000., 4) for t in df_curr.index]
        ax.bar(t, df_curr['Time'], align='center', width=0.05)
#        ax.set_xlabel('$s$')
#        ax.set_ylabel('Spikes')
    return fig, axes


def recordElectrodeEnviroment(df, posX, posY, dX, dY):
    """
    Shows firing Rate for Neurons at given Position (posX, posY) and its enivroment (posX+dX, posY+dY)
    :return: Plot
    """
    times = [t for (x, y), t in zip(df['Position'], df['Time']) if posX-dX <= x <= posX+dX and posY-dY <= y <= posY+dY]
    figure = plt.hist(times, bins=50)
    return figure


class RandomBalancedNetwork:
    def __init__(self, parameters):
        self.parameters = parameters
        self.gridSize = parameters['Columns']*parameters['Rows']
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": 0.1, "print_time": True,
                              "overwrite_files": True,
                              "local_num_threads": 6
                             })
        nest.CopyModel('iaf_psc_alpha', 'exci')
        nest.CopyModel('iaf_psc_alpha', 'inhi')
        nest.CopyModel('static_synapse', 'exc', {'weight': self.parameters['Excitational Weight']})
        nest.CopyModel('static_synapse', 'inh', {'weight': self.parameters['Inhibitory Weight']})
        # nest.CopyModel('static_synapse', 'exc')
        # nest.CopyModel('static_synapse', 'inh')
        nest.CopyModel('static_synapse', 'inh_strong', {'weight': self.parameters['Weight Stimulus']})
        self.l = tp.CreateLayer({'rows': self.parameters['Rows'],
                                 'columns': self.parameters['Columns'],
                                 'elements': ['exci', self.parameters['Number excitational cells'],
                                              'inhi', self.parameters['Number inhibitory cells']],
                                 'edge_wrap': False})
        cdict_e2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius excitational']}},
                     'kernel': {'gaussian': {'p_center': parameters['Excitatory Weight Maximum'], 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': parameters['e2i delay'], 'a': parameters['e2i delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'exc',
                     'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
        cdict_e2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius excitational']}},
                     'kernel': {'gaussian': {'p_center': parameters['Excitatory Weight Maximum'], 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': parameters['e2e delay'], 'a': parameters['e2e delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'exc',
                     'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
        cdict_i2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': parameters['Inhibitory Weight Maximum'], 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': parameters['i2e delay'], 'a': parameters['i2e delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'inh',
                     'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
        cdict_i2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': parameters['Inhibitory Weight Maximum'], 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': parameters['i2i delay'], 'a': parameters['i2i delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'inh',
                     'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
        tp.ConnectLayers(self.l, self.l, cdict_e2i)
        tp.ConnectLayers(self.l, self.l, cdict_e2e)
        tp.ConnectLayers(self.l, self.l, cdict_i2i)
        tp.ConnectLayers(self.l, self.l, cdict_i2e)
        stim = tp.CreateLayer({'rows': 1,
                               'columns': 1,
                               'elements': 'poisson_generator'})
        stim_i = nest.GetLeaves(stim, local_only=True)[0]
        nest.SetStatus(stim_i, {'rate': parameters['Background rate']})
        background_stim_dict = {'connection_type': 'divergent',
                                # 'mask': {'grid': {'rows': self.parameters['Rows'],
                                                  # 'columns': self.parameters['Columns']}},
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
        self.df_ex = makePandas(self.events_ex, tp.FindCenterElement(self.l)[0])
        self.df_in = makePandas(self.events_in, tp.FindCenterElement(self.l)[0])

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

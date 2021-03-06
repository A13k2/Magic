import nest
# Set verbosity to send warnings or higher
# M_ALL=0,  display all messages
# M_DEBUG=5,  display debugging messages and above
# M_STATUS=7,  display status messages and above
# M_INFO=10, display information messages and above
# M_DEPRECATED=18, display deprecation warnings and above
# M_WARNING=20, display warning messages and above
# M_ERROR=30, display error messages and above
# M_FATAL=40, display failure messages and above
# M_QUIET=100, suppress all messages
nest.set_verbosity(20)
import os
import nest.topology as tp
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import topology_2d_helper as magic
from progressbar import *


"""
------------------------------------------
Simulation Control, what to start when....
------------------------------------------
"""
def simulationAndAnalysis(parameters, curr_folder='.'):
    delay_visualisation_linear(parameters, curr_folder+'/delay.png')
    weightVisualisation(parameters, curr_folder+'/weights.pdf')
    simulation = magic.RandomBalancedNetwork(parameters)
    simulation.start_simulation()
    simulation.writeParametersToFile(curr_folder + '/parameters.txt')
    simulation.pickle_dump(curr_folder)
    # makeDir(curr_folder+'/excitatory_neurons')
    # makeDir(curr_folder+'/inhibitory_neurons')
    # stimulationControlLazy(simulation, curr_folder)
    # #clusteringPlotLazy(simulation, parameters, curr_folder)
    # #fanoFactorTimeLazy(simulation, curr_folder)
    # #recordElectrodeEnviromentLazy(simulation, curr_folder)
    # recordElectrodeLazy(simulation, curr_folder)
    # spikeCountHistogramLazy(simulation, curr_folder)
    # rasterPlotLazy(simulation, curr_folder)
    # distancePlotsLazy(1000., 2000., 500., simulation, curr_folder)

def simulate_and_dump(parameters, folder):
    magic.makeDir(folder)
    simulation = magic.RandomBalancedNetwork(parameters)
    simulation.start_simulation()
    simulation.writeParametersToFile(folder + '/parameters.txt')
    simulation.pickle_dump(folder)


def average_firng_rates(parameters, curr_folder='.'):
    """
    simulates network with parameters and saves values in pickle file
    """
    import pickle
    parameters['Time of stimulation'] = 1000.
    num_i_neurons, num_e_neurons = num_of_neurons(parameters)
    for target in ['exci', 'inhi']:
        for jee in [1.]:
            for back_i in [0., 300., 600., 1200.]:
                for stim_weight in [-300000., 300000.]:
                    parameters['target'] = target
                    parameters['Jee'] = jee
                    parameters['Weight Stimulus'] = stim_weight
                    parameters['Background rate inhibitory'] = back_i
                    # simulation = magic.RandomBalancedNetwork(parameters)
                    print(target, jee, back_i, stim_weight)
                    print('%s/stimulus/%s_%d_%02d_%d.p' % ( curr_folder,
                        target[:-1], jee, int(back_i/100.), stim_weight))
                    # simulation.start_simulation()
                    # average_e = firing_rate_time(simulation.df_ex, num_e_neurons)
                    # average_i = firing_rate_time(simulation.df_in, num_i_neurons)
                    # with open('%s/stimulus/%s_%d_%02d_%d.p' % (
                        # curr_folder, target[:-1], jee, int(back_i/100.), stim_weight
                        # ), 'wb') as f:
                    # # with open(curr_folder + '/stimulus/' + target[:-1] + '_' +
                            # # str(int(jee))+'_'+str(int(back_i/100.))+'_'+str(int(stim_weight))+'.p','wb') as f:
                        # pickle.dump((average_e, average_i), f)

def firing_rate_time(df, num_neurons, num_bins=20):
    bins = np.linspace(df['Time'].min(), df['Time'].max(), num_bins)
    groups = df.groupby(pd.cut(df['Time'], bins))
    groups = groups['Time'].count()/num_neurons
    for g, ind in zip(groups, groups.index):
        groups[ind] = g/(ind.right-ind.left)*1000.
    return groups

def num_of_neurons(parameters):
    grid_size = parameters['Columns']*parameters['Rows']
    return (parameters['Number inhibitory cells']*grid_size, parameters['Number excitational cells']*grid_size)

def tsodyks_analysis(parameters, curr_folder='.'):
    """
    Plot average firing rate of populations E and I for different noise rate of I
    Saves different Plots for different J_ee
    """
    def e_i_of_i_ext():
        i_start = 0.
        i_max = 30000.
        i_steps = 10
        i_ext = np.linspace(i_start, i_max, i_steps)
        i_average = []
        e_average = []
        grid_points = parameters['Rows']*parameters['Columns']
        num_i_neurons = parameters['Number inhibitory cells']*grid_points
        num_e_neurons = parameters['Number excitational cells']*grid_points
        t_min = parameters['Time before stimulation']+parameters['Time of stimulation']
        t_max = t_min+parameters['Time after Stimulation']
        for inh_background in i_ext:
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('Running "Tsodyks analysis" with currently i_{ext} = '
                  + str(inh_background) + ' of i_{ext, max} = ' + str(i_max))
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            parameters['Background rate inhibitory'] = inh_background
            simulation = magic.RandomBalancedNetwork(parameters)
            simulation.start_simulation()
            simulation.writeParametersToFile(curr_folder+'/parameters.txt')
            i_average.append(average_firing_rate(simulation.df_in, num_i_neurons, t_min, t_max))
            e_average.append(average_firing_rate(simulation.df_ex, num_e_neurons, t_min, t_max))
            print('Average firing rate for inhibitory population: '+
                  str(i_average[-1]))
            print('Average firing rate for excitatory population: '+
                  str(e_average[-1]))
        return i_average, e_average, i_ext
    delay_visualisation_linear(parameters, curr_folder+'/delay.pdf')
    weightVisualisation(parameters, curr_folder+'/weights.pdf')
    j_ee_min = 3
    j_ee_num = 24
    j_ee_max = 12
    j_ee = np.linspace(j_ee_min, j_ee_max, j_ee_num)
    for j_ in j_ee:
        # parameters['Jei'] = j_
        parameters['Jee'] = j_
        i_average, e_average, i_ext = e_i_of_i_ext()
        plt.clf()
        plt.plot(i_ext, i_average, label=r'$\hat{I}$')
        plt.plot(i_ext, e_average, label=r'$\hat{E}$')
        plt.xlabel(r'$i_{ext}$')
        plt.ylabel(r'$\hat{\nu}$')
        plt.legend()
        plt.savefig(curr_folder+'/IE_vs_i_ext_j_ee_'+str(round(j_,1))+'.png')

def tsodyks_analysis_phase_plane(parameters, curr_folder='.', e_range=np.arange(0., 300., 15.),
            i_range=np.arange(0., 300., 15.), j_ee_range=np.arange(1.,9.,2.)):
    """
    Plot quiver Plot for different noise of E and I
    """
    import pickle
    sim_size = 0
    it = 0
    pbar = ProgressBar()
    grid_points = parameters['Rows']*parameters['Columns']
    num_i_neurons = parameters['Number inhibitory cells']*grid_points
    num_e_neurons = parameters['Number excitational cells']*grid_points
    t_min = parameters['Time before stimulation']+parameters['Time of stimulation']
    t_max = t_min+parameters['Time after Stimulation']

    @np.vectorize
    def e_i(e_ext, i_ext, paramters):
        """
        Calculate average firing rate of population E and I for given external
        noise
        """
        parameters['Background rate inhibitory'] = i_ext
        parameters['Background rate excitatory'] = e_ext
        simulation = magic.RandomBalancedNetwork(parameters)
        simulation.start_simulation()
        e_average = average_firing_rate(simulation.df_ex, num_e_neurons, t_min, t_max)
        i_average = average_firing_rate(simulation.df_in, num_i_neurons, t_min, t_max)
        nonlocal it
        pbar.update(it)
        it += 1
        return e_average, i_average

    def calculate_from_mesh(e_arange, i_arange):
        """
        Calcualtes mesh grid for e_arange and i_arange
        calculates return value for each grid point
        returns E,I,e_phase_plane_analysis/new/average,i_average
        """
        E, I = np.meshgrid(e_arange, i_arange)
        nonlocal sim_size
        sim_size = E.size
        pickle.dump((E,I), open(curr_folder+'E_I.p', 'wb'))
        pickle.dump(parameters, open(curr_folder+'parameter.p', 'wb'))
        for j_ee in j_ee_range:
            print("")
            print("Starting simulations for J_ee = " + str(j_ee))
            nonlocal pbar
            pbar = ProgressBar(widgets=[Percentage(), ' ', ETA(), ' ', Bar()], maxval=sim_size).start()
            parameters['Jee'] = j_ee
            E_average, I_average = e_i(E, I, parameters)
            pickle.dump((E_average, I_average), open(curr_folder+'e_i_jee_'+str(j_ee).replace('.','_')+'.p', 'wb'))
            nonlocal it
            it = 0

    calculate_from_mesh(e_range, i_range)


"""
--------------
General helper
--------------
"""
def average_firing_rate(df, num_neurons, t_min, t_max):
    df = df[df['Time'] >= t_min]
    df = df[df['Time'] <= t_max]
    return len(df)/float(num_neurons*((t_max-t_min)/1000.))


def makeDir(folder):
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


def makePandas(senderEvents, ctr):
    """
    Return a pandas DataFrame (Sender, Time, Position, Distance from ctr) for the recordings of a spike detector
    :param senderEvents: Recordings from spike detector
    :param ctr: Center Element(GID)
    :return: Pandas dataframe with Sender, Time, Position and Distance
    """
    senders = [int(i) for i in senderEvents['senders']]  #int necessary, because GID's for tp.GetPosition. Otherwise error
    times = [i for i in senderEvents['times']]
    positions = [i for i in tp.GetPosition(senders)]
    if len(senders) >= 2:
        distances = [i for i in tp.Distance(senders, [ctr])]
        return pd.DataFrame({'Sender': senders,
                             'Time': times, 'Distance from Center': distances,
                             'Position': positions
                             })
    else:
        return pd.DataFrame({'Sender': senders,
                             'Time': times,
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


"""
-----------
Math helper
-----------
"""
def linear(x,a,b):
    return a+x*b


def gaussian(x, mu, sig, maximum=1.):
    return maximum*np.exp(-np.power(x - mu, 2.) / (2. * np.power(sig, 2.)))


"""
-----------
Plot Helper
-----------
"""
def networkVisualizationLazy(parameters, folder):
    '''
    Plots topology, weights and more in a subfolder called network_visu
    '''
    delay_visualisation_linear(parameters, folder+'/delay.pdf')
    weightVisualisation(parameters, folder+'/weights.pdf')
    simulation = magic.RandomBalancedNetwork(parameters)
    simulation.plotKernel(folder)


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


"""
--------------
Visualisations
--------------
"""
def weightVisualisation(parameters, file):
    plt.clf()
    x = np.arange(-0.5,0.51,0.01)
    plt.plot(x, gaussian(x, 0., parameters['Sigma excitational'], maximum=parameters['Jee Connectivity']), label='Excitational')
    plt.plot(x, gaussian(x, 0., parameters['Sigma inhibitory'], maximum=parameters['Jii Connectivity']), label='Inhibitory')
    plt.plot(x, gaussian(x, 0., parameters['Sigma excitational'], maximum=parameters['Jee Connectivity'])-gaussian(x, 0., parameters['Sigma inhibitory'], maximum=parameters['Jii Connectivity']), label='Excitational - Inhibitory')
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

def distance_firing_rate(df, number_of_bins, title):
    duration = df['Time'].values.max() - df['Time'].values.min()
    print(duration)
    grouped = df.groupby([pd.cut(df['Distance from Center'], number_of_bins), pd.cut(df['Time'], number_of_bins)])
    neurons_per_area = df.groupby(pd.cut(df['Distance from Center'], number_of_bins))['Sender'].nunique()
#     Magic to change index from categorical to right side of interval
    tmp = grouped.count().unstack().divide(neurons_per_area, axis=0).multiply(1000./(duration/number_of_bins))['Distance from Center']
    tmp.index = [a.right for a in tmp.index]
    ax = tmp.plot(title=title)
    ax.set_ylabel(r'$\hat{\nu}$')
    return ax

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
    x, y = spike_count_histogram(df, t_start, t_stop, t_step, number_of_neurons, grid_size)
    return plt.bar(x[:-1], y, width=0.9*t_step)



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
        tauSyn = 2.5  # synaptic time constant in ms
        tauMem = 20.0 # time constant of membrane potential in ms
        CMem = 250.0  # capacitance of membrane in in pF
        theta = 20.0  # membrane threshold potential in mV
        neuron_params = {"C_m": CMem,
                 "tau_m": tauMem,
                 "tau_syn_ex": tauSyn,
                 "tau_syn_in": tauSyn,
                 "t_ref": 2.0,
                 "E_L": 0.0,
                 "V_reset": 0.0,
                 "V_m": 0.0,
                 "V_th": theta}
        self.parameters = parameters
        self.gridSize = parameters['Columns']*parameters['Rows']
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": 0.1, "print_time": False,
                              "overwrite_files": True,
                              "local_num_threads": 6
                             })
        nest.CopyModel('iaf_psc_alpha', 'exci')
        nest.SetDefaults('exci', neuron_params)
        nest.CopyModel('iaf_psc_alpha', 'inhi')
        nest.SetDefaults('inhi', neuron_params)
        nest.CopyModel('static_synapse', 'exc', {'weight': self.parameters['Jei']})
        nest.CopyModel('static_synapse', 'background', {'weight': self.parameters['Background weight']})
        nest.CopyModel('static_synapse', 'inh', {'weight': self.parameters['Jie']})
        nest.CopyModel('static_synapse', 'exc_recurrent', {'weight': self.parameters['Jee']})
        nest.CopyModel('static_synapse', 'inh_recurrent', {'weight': self.parameters['Jii']})
        nest.CopyModel('static_synapse', 'inh_strong', {'weight': self.parameters['Weight Stimulus']})
        self.l = tp.CreateLayer({'rows': self.parameters['Rows'],
                                 'columns': self.parameters['Columns'],
                                 'elements': ['exci', self.parameters['Number excitational cells'],
                                              'inhi', self.parameters['Number inhibitory cells']],
                                 'edge_wrap': False})
        self.cdict_e2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius excitational']}},
                     'kernel': {'gaussian': {'p_center': parameters['Jei Connectivity'], 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': parameters['e2i delay'], 'a': parameters['e2i delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'exc'}
                     #'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
                     #'synapse_model': 'exc',
                     #'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
        self.cdict_e2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius excitational']}},
                     'kernel': {'gaussian': {'p_center': parameters['Jee Connectivity'], 'sigma': self.parameters['Sigma excitational']}},
                     'delays': {'linear': {'c': parameters['e2e delay'], 'a': parameters['e2e delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'exc_recurrent'}
                     #'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
                     #'synapse_model': 'exc_recurrent',
                     #'weights': {'uniform': {'min': parameters['Excitational Weight']*0.2, 'max': parameters['Excitational Weight']}}}
        self.cdict_i2e = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': parameters['Jie Connectivity'], 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': parameters['i2e delay'], 'a': parameters['i2e delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'exci'},
                     'synapse_model': 'inh'}
                     #'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
                     #'synapse_model': 'inh',
                     #'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
        self.cdict_i2i = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': self.parameters['Radius inhibitory']}},
                     'kernel': {'gaussian': {'p_center': parameters['Jii Connectivity'], 'sigma': self.parameters['Sigma inhibitory']}},
                     'delays': {'linear': {'c': parameters['i2i delay'], 'a': parameters['i2i delay']*parameters['delay growth multiplier']}},
                     'sources': {'model': 'inhi'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'inh_recurrent'}
                     # 'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
                     # 'synapse_model': 'inh_recurrent',
                     # 'weights': {'uniform': {'max': parameters['Inhibitory Weight']*0.2, 'min': parameters['Inhibitory Weight']}}}
        tp.ConnectLayers(self.l, self.l, self.cdict_e2i)
        tp.ConnectLayers(self.l, self.l, self.cdict_e2e)
        tp.ConnectLayers(self.l, self.l, self.cdict_i2i)
        tp.ConnectLayers(self.l, self.l, self.cdict_i2e)
        """
        Creating Background
        """
        #Excitatory Background
        background_e = tp.CreateLayer({'rows': 1,
                                       'columns': 1,
                                       'elements': 'poisson_generator'})
        background_e_leaves = nest.GetLeaves(background_e, local_only=True)[0]
        nest.SetStatus(background_e_leaves, {'rate': parameters['Background rate excitatory']})
        background_e_dict = {'connection_type': 'divergent',
                             'targets': {'model': 'exci'},
                             'synapse_model': 'background'}
        tp.ConnectLayers(background_e, self.l, background_e_dict)
        # Inhibitory Background
        background_i = tp.CreateLayer({'rows': 1,
                                       'columns': 1,
                                       'elements': 'poisson_generator'})
        background_i_leaves = nest.GetLeaves(background_i, local_only=True)[0]
        nest.SetStatus(background_i_leaves, {'rate': parameters['Background rate inhibitory']})
        background_i_dict = {'connection_type': 'divergent',
                             'targets': {'model': 'inhi'},
                             'synapse_model': 'background'}
        tp.ConnectLayers(background_i, self.l, background_i_dict)
        # Excitation/ Inhibition in center
        stim2 = tp.CreateLayer({'rows': 1,
                                'columns': 1,
                                'elements': 'poisson_generator'})
        self.stim2_i = nest.GetLeaves(stim2, local_only=True)[0]
        nest.SetStatus(self.stim2_i, {'rate': 0.0})
        self.cdict_stim2 = {'connection_type': 'divergent',
                            'kernel': {'gaussian': {'p_center': 1., 'sigma': self.parameters['Sigma Stimulus']}},
                            'mask': {'circular': {'radius': self.parameters['Radius stimulus']},
                                     'anchor': [0., 0.]},
                            # 'targets': {'model': 'inhi'},  # set population
                            # to be affected by stimulus here
                            # 'targets': {'model': 'exci'},
                            'targets': {'model': self.parameters['target']},
                            'synapse_model': 'inh_strong'}
        tp.ConnectLayers(stim2, self.l, self.cdict_stim2)
        # Connect only one neurons per node to spike detector
        self.rec_ex = nest.Create("spike_detector")
        self.rec_in = nest.Create("spike_detector")
        grid_points = self.parameters['Columns']*self.parameters['Rows']
        number_exc_cells = grid_points*self.parameters['Number excitational cells']
        # layer_exc = tuple(i for i in range(2, grid_points+2))
        # layer_inh = tuple(i for i in range(number_exc_cells+2, number_exc_cells+grid_points+2))
        # nest.Connect(layer_exc, self.rec_ex)
        # nest.Connect(layer_inh, self.rec_in)
        # Connect all neurons per node to spike detector
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

    def plotKernel(self, folder):
        fig, ax = plt.subplots()
        tp.PlotLayer(self.l, fig, nodesize=10)
        ctr_elem = tp.FindCenterElement(self.l)
        # import ipdb; ipdb.set_trace()
        # tp.PlotKernel(ax, ctr_elem, mask=self.cdict_i2i['mask'], kern=self.cdict_i2i['kernel'])
        tp.PlotTargets(ctr_elem, self.l, fig=fig, tgt_color='red',
                mask=self.cdict_e2e['mask'], kernel=self.cdict_e2e['kernel'],
                kernel_color='green', mask_color='green',
                syn_type='exc_recurrent')
        # tp.PlotKernel(ax, ctr_elem, mask=self.cdict_e2e['mask'], mask_color='green')
        plt.xlabel(r'X')
        plt.ylabel(r'Y')
        plt.xlim(-.5, .5)
        plt.ylim(-.5, .5)
        plt.savefig(folder+'/kernel.pdf')

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
        # self.events_ex = nest.GetStatus(self.rec_ex, "events")[0]
        # self.events_in = nest.GetStatus(self.rec_in, "events")[0]
        self.df_ex = makePandas(self.events_ex, tp.FindCenterElement(self.l)[0])
        self.df_in = makePandas(self.events_in, tp.FindCenterElement(self.l)[0])

    def pickle_dump(self, folder):
        import pickle
        with open(folder+'/parameter.pickle', 'wb') as f:
            pickle.dump(self.parameters, f)
        self.df_ex.to_pickle(folder+'/df_ex.pickle')
        self.df_in.to_pickle(folder+'/df_in.pickle')
        with open(folder+'/rec.pickle','wb') as f:
            pickle.dump((self.events_ex, self.events_in), f)

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

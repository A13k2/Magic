# import for missing pythonpath
# import imp
# nest= imp.load_source('NEST', '/cm/shared/software/NEST/2.14.0-foss-2016b-Python-3.6.1/lib64/python3.6/site-packages/nest/__init__.py')

import general_helper as gen
import nest.topology as tp
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import topology_2d_helper as magic

base_folder = '/home/adrossel/Magic/figures/'
date_folder = '18_08_01/'
sub_folder = 'test/'
all_folder = base_folder + date_folder + sub_folder
pd.set_option('display.max_columns', 500)


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
        magic.visualization(df_before_stim, title + ' Neurons before stimulation, ' + network.parameters['Name'])
        plt.savefig(folder + neurons_folder + '/visu_before_stim_' + title + '_neurons.png',
                    dpi=300)
        df_while_stim = df[df['Time'] > network.parameters['Time before stimulation']]
        df_while_stim = df_while_stim[
            df_while_stim['Time'] <= network.parameters['Time before stimulation'] + network.parameters[
                'Time of stimulation']]
        plt.close()
        magic.visualization(df_while_stim, title + ' Neurons while stimulation, ' + network.parameters['Name'])
        plt.savefig(folder + neurons_folder + '/visu_while_stim_' + title + '_neurons.png',
                    dpi=300)
        df_after_stim = df[
            df['Time'] > network.parameters['Time before stimulation'] + network.parameters['Time of stimulation']]
        plt.close()
        magic.visualization(df_after_stim, title + ' Neurons after stimulation, ' + network.parameters['Name'])
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


def automation():
    parameters = {'Name': 'inhibition',
                  'Columns': 40,
                  'Rows': 40,
                  'Excitational Weight': 5.0,
                  'Radius excitational': 0.1,
                  'Sigma excitational': 0.05,
                  'Inhibitory Weight': -5.0,
                  'Radius inhibitory': 0.15,
                  'Sigma inhibitory': 0.075,
                  'Number excitational cells': 5,
                  'Number inhibitory cells': 5,
                  'Weight Stimulus': -3000.,
                  'Radius stimulus': 0.1,
                  'Sigma Stimulus': 0.05,
                  'Stimulus rate': 40000.,
                  'Background rate': 35000.,
                  'Time before stimulation': 1000.,
                  'Time of stimulation': 500.,
                  'Time after Stimulation': 1000.,
                  }
    radius_inhibs = [0.1, 0.2, 0.2, 0.1, 0.2]
    sigma_inhibs = [0.075, 0.1, 0.1, 0.05, 0.75]
    radius_excis = [0.1, 0.2, 0.1, 0.2, 0.1]
    sigma_excis = [0.05, 0.1, 0.05, 0.1, 0.5]
    weight_inhis = [-5., -5.0, -4.0, -10.0, -20.0]
    weight_excis = [5., 5.0, 4.0, 10.0, 20.0]
    sim_folder = base_folder + date_folder
    for cols_rows in [60, 80]:
        col_folder = sim_folder + '/colsRows_' + str(cols_rows)
        gen.create_folder(col_folder)
        parameters['Columns'] = cols_rows
        parameters['Rows'] = cols_rows
        for background_rate in [15000., 20000., 25000., 30000., 35000., 40000.]:
            background_folder = col_folder + '/background_' + str(background_rate)
            gen.create_folder(background_folder)
            parameters['Background rate'] = background_rate
            for radius_inhib, sigma_inhib, radius_exci, sigma_exci, weight_inhi, weight_exci  in zip(radius_inhibs, sigma_inhibs, radius_excis, sigma_excis, weight_inhis, weight_excis):
                name = "inh_(" + str(radius_inhib) + "," + str(sigma_inhib) + ")_exci_(" + str(radius_exci) + "," + str(sigma_exci) + "_wight_inhi_" + str(weight_inhi) + "weight_exci_" + str(weight_exci) + ")"
                curr_folder = background_folder + '/' + name
                if not os.path.exists(curr_folder):
                    os.mkdir(curr_folder)
                    os.mkdir(curr_folder+'/excitatory_neurons')
                    os.mkdir(curr_folder+'/inhibitory_neurons')
                parameters['Name'] = name
                parameters['Radius excitational'] = radius_exci
                parameters['Radius inhibitory'] = radius_inhib
                parameters['Sigma excitational'] = sigma_exci
                parameters['Sigma inhibitory'] = sigma_inhib
                simulation = magic.RandomBalancedNetwork(parameters)
#                simulations.append(simulation)
                simulation.start_simulation()
                simulation.writeParametersToFile(curr_folder + '/parameters.txt')
#                recordElectrodeLazy(simulation, curr_folder)
                recordElectrodeEnviromentLazy(simulation, curr_folder)
                fanoFactorTimeLazy(simulation, curr_folder)
                spikeCountHistogramLazy(simulation, curr_folder)
                stimulationControlLazy(simulation, curr_folder)
                rasterPlotLazy(simulation, curr_folder)
#                distancePlotsLazy(0., 0.5, 0.05, simulation, curr_folder)
                del simulation


def testSimulation():
    parameters = {'Name': 'inhibition',
                  'Columns': 40,
                  'Rows': 40,
                  'Excitational Weight': 5.0,
                  'Radius excitational': 0.1,
                  'Sigma excitational': 0.05,
                  'Inhibitory Weight': -5.0,
                  'Radius inhibitory': 0.15,
                  'Sigma inhibitory': 0.075,
                  'Number excitational cells': 5,
                  'Number inhibitory cells': 5,
                  'Weight Stimulus': -3000.,
                  'Radius stimulus': 0.1,
                  'Sigma Stimulus': 0.05,
                  'Stimulus rate': 40000.,
                  'Background rate': 35000.,
                  'Time before stimulation': 1000.,
                  'Time of stimulation': 500.,
                  'Time after Stimulation': 1000.,
                  }
    simulation = magic.RandomBalancedNetwork(parameters)
    simulation.start_simulation()
    simulation.writeParametersToFile(all_folder + 'parameters.txt')
    neurons_x_position = tp.FindNearestElement(simulation.l, (0.,0.))
    neurons_position = tp.GetPosition(neurons_x_position)
    magic.recordElectrodeEnviroment(simulation.df_ex, neurons_position[0][0], neurons_position[0][1], 0.1, 0.1)
    plt.savefig(all_folder + '/excitatory_neurons/electrode_' + str(round(neurons_position[0], 4)) + ', ' + str(round(neurons_position[1], 4)) + '.png', dpi=300)
    plt.close()


automation()
#testSimulation()

#parameters = {'Name': 'inhibition',
#              'Columns': 40,
#              'Rows': 40,
#              'Excitational Weight': 5.0,
#              'Radius excitational': 0.1,
#              'Sigma excitational': 0.05,
#              'Inhibitory Weight': -5.0,
#              'Radius inhibitory': 0.15,
#              'Sigma inhibitory': 0.075,
#              'Number excitational cells': 5,
#              'Number inhibitory cells': 5,
#              'Weight Stimulus': -3000.,
#              'Radius stimulus': 0.1,
#              'Sigma Stimulus': 0.05,
#              'Stimulus rate': 40000.,
#              'Background rate': 35000.,
#              'Time before stimulation': 1000.,
#              'Time of stimulation': 500.,
#              'Time after Stimulation': 1000.,
#              }
#
#
#exc = magic.RandomBalancedNetwork(parameters)
#exc.start_simulation()
#
#
#exc.writeParametersToFile(all_folder + 'parameters.txt')
#recordElectrodeLazy(exc)
#fanoFactorTimeLazy(exc)
#spikeCountHistogramLazy(exc)
#stimulationControlLazy(exc)
## distancePlotsLazy(0.,.5,.1,exc)
#magic.visualization(exc.df_ex, "Test")
#
#magic.distanceFiringRate(exc.events_ex, tp.FindCenterElement(exc.l)[0],
#                         neurons_per_gridpoint=exc.parameters['Number inhibitory cells'], title="Excitatory",
#                         gridSize=exc.gridSize)
#plt.savefig(all_folder + 'dist_exci_exci.png', dpi=300)
#magic.distanceFiringRate(exc.events_in, tp.FindCenterElement(exc.l)[0],
#                         neurons_per_gridpoint=exc.parameters['Number excitational cells'], title='Inhibitory',
#                         gridSize=exc.gridSize)
#plt.savefig(all_folder + 'dist_inhi_exci.png', dpi=300)
#
## fig = tp.PlotLayer(l)
## ctr = tp.FindCenterElement(l)
## tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
## plt.show()
## fig = tp.PlotLayer(l)
## ctr = tp.FindCenterElement(l)
## tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')
#
#
# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
# plt.show()
# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')
# plt.show()

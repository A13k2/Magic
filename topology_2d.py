import nest
import nest.topology as tp
import numpy as np
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import matplotlib.animation as animation
from pandas.util.testing import network

import topology_2d_helper as magic

base_folder = '/home/alex/Magic/figures/'
date_folder = '18_07_06/'
sub_folder = 'inhibition40x40/'
all_folder = base_folder + date_folder + sub_folder

def distancePlotsLazy(Start, End, Step, network):
    for events, title in zip([network.events_ex, network.events_in], ["Excitatory", "Inhibitory"]):
        for start, end in zip(np.arange(Start, End, Step), np.arange(Start + Step, End + Step, Step)):
            plt.clf()
            currentTitle = title + 'Neurons from ' + str(start) + ' to ' + str(end)
            print(currentTitle)
            times, neurons = magic.distance(events, start, end, tp.FindCenterElement(network.l)[0])
            magic.raster_plot(senders=times, timeS=neurons, title=currentTitle, gridSize=network.gridSize)
            plt.savefig(
                all_folder + 'dist_' + str(end) + '_raster_' + title + '_' + network.parameters['Name'] + '.png',
                dpi=300)


def stimulationControlLazy(network):
    for df, title in zip([network.df_ex, network.df_in], ['Excitatory', 'Inhibitory']):
        df_before_stim = df[df['Time'] <= network.parameters['Time before stimulation']]
        plt.clf()
        magic.visualization(df_before_stim, title + ' Neurons before stimulation, ' + network.parameters['Name'])
        plt.savefig(all_folder + 'visu_before_stim_' + title + '_neurons_' + network.parameters['Name'] + '.png',
                    dpi=300)
        df_while_stim = df[df['Time'] > network.parameters['Time before stimulation']]
        df_while_stim = df_while_stim[
            df_while_stim['Time'] <= network.parameters['Time before stimulation'] + network.parameters[
                'Time of stimulation']]
        plt.clf()
        magic.visualization(df_while_stim, title + ' Neurons while stimulation, ' + network.parameters['Name'])
        plt.savefig(all_folder + 'visu_while_stim_' + title + '_neurons_' + network.parameters['Name'] + '.png',
                    dpi=300)
        df_after_stim = df[
            df['Time'] > network.parameters['Time before stimulation'] + network.parameters['Time of stimulation']]
        plt.clf()
        magic.visualization(df_after_stim, title + ' Neurons after stimulation, ' + network.parameters['Name'])
        plt.savefig(all_folder + 'visu_after_stim_' + title + '_neurons_' + network.parameters['Name'] + '.png',
                    dpi=300)


def fanoFactorTimeLazy(network):
    plt.clf()
    magic.fanoFactorTimePlot(network.df_ex, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number excitational cells'], network.gridSize, bins=20)
    plt.savefig(all_folder + 'fano_factor_time_excitatory_neurons_' + network.parameters['Name'] + '.png',
                dpi=300)
    plt.clf()
    magic.fanoFactorTimePlot(network.df_in, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number inhibitory cells'], network.gridSize, bins=20)
    plt.savefig(all_folder + 'fano_factor_time_inhibitory_neurons_' + network.parameters['Name'] + '.png',
                dpi=300)


def spikeCountHistogramLazy(network):
    plt.clf()
    magic.spike_count_histogram_plot(network.df_ex, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number excitational cells'], network.gridSize)
    plt.savefig(all_folder + 'spike_count_histogram_excitatory_neurons_' + network.parameters['Name'] + '.png',
                dpi=300)
    plt.clf()
    magic.spike_count_histogram_plot(network.df_in, 0., (network.parameters['Time before stimulation']+network.parameters['Time of stimulation']+network.parameters['Time after Stimulation'])/1000., .01, network.parameters['Number inhibitory cells'], network.gridSize)
    plt.savefig(all_folder + 'spike_count_histogram_inhibitory_neurons_' + network.parameters['Name'] + '.png',
                dpi=300)


def recordElectrodeLazy(network):
    neurons_x_position_distance = 0.5/(network.parameters['Columns']/2.)
    neurons_x_position = [(i, 0.) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)]
    neurons_x_position += ([(0., i) for i in np.arange(0., 0.5, neurons_x_position_distance*2.)])
    neurons_x_position = tp.FindNearestElement(network.l, neurons_x_position)
    neurons_position = tp.GetPosition(neurons_x_position)
    for position in neurons_position:
        plt.clf()
        magic.recordElectrode(exc.df_ex, position[0], position[1], 8)
        plt.savefig(all_folder + 'record_electrode_exci_neurons' + str(round(position[0],4)) + ', ' + str(round(position[1], 4)) + network.parameters['Name'] + '.png',
                    dpi=300)
        plt.clf()
        magic.recordElectrode(exc.df_in, position[0], position[1], 8)
        plt.savefig(all_folder + 'record_electrode_inhi_neurons_' + str(round(position[0],4)) + ', ' + str(round(position[1], 4)) + network.parameters['Name'] + '.png',
                    dpi=300)


parameters = {'Name': 'inhibition',
              'Columns': 40,
              'Rows': 40,
              'Excitational Weight': 5.0,
              'Radius excitational': 0.35,
              'Sigma excitational': 0.15,
              'Inhibitory Weight': -5.0,
              'Radius inhibitory': 0.5,
              'Sigma inhibitory': 0.18,
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

pd.set_option('display.max_columns', 500)

exc = magic.RandomBalancedNetwork(parameters)
exc.start_simulation()


exc.writeParametersToFile(all_folder + 'parameters.txt')
recordElectrodeLazy(exc)
fanoFactorTimeLazy(exc)
spikeCountHistogramLazy(exc)
stimulationControlLazy(exc)
# distancePlotsLazy(0.,.5,.1,exc)
magic.visualization(exc.df_ex, "Test")
magic.raster_plot(eventSenders=exc.events_ex, gridSize=exc.gridSize)
plt.savefig(all_folder + 'raster_exci_' + exc.parameters['Name'] + '.png', dpi=300)
magic.raster_plot(eventSenders=exc.events_in, gridSize=exc.gridSize)
plt.savefig(all_folder + 'raster_inhi_' + exc.parameters['Name'] + '.png', dpi=300)

magic.distanceFiringRate(exc.events_ex, tp.FindCenterElement(exc.l)[0],
                         neurons_per_gridpoint=exc.parameters['Number inhibitory cells'], title="Excitatory",
                         gridSize=exc.gridSize)
plt.savefig(all_folder + 'dist_exci_exci.png', dpi=300)
magic.distanceFiringRate(exc.events_in, tp.FindCenterElement(exc.l)[0],
                         neurons_per_gridpoint=exc.parameters['Number excitational cells'], title='Inhibitory',
                         gridSize=exc.gridSize)
plt.savefig(all_folder + 'dist_inhi_exci.png', dpi=300)

# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
# plt.show()
# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')


# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
# plt.show()
# fig = tp.PlotLayer(l)
# ctr = tp.FindCenterElement(l)
# tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')
# plt.show()

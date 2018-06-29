import nest
import nest.topology as tp
import numpy as np
import matplotlib.pyplot as plt
import nest.raster_plot
import pandas as pd
import matplotlib.animation as animation
import topology_2d_helper as magic




cols_rows=80
gridSize = cols_rows**2
base_folder = '/home/alex/compu/figures/'
date_folder = '18_06_29/'
sub_folder = 'test_write_params/'
all_folder = base_folder+date_folder+sub_folder

parameters = {'Columns': 80,
              'Rows': 80,
              'Excitational Weight': 5.0,
              'Radius excitational': 0.05,
              'Sigma excitational': 0.025,
              'Inhibitory Weight': -20.0,
              'Radius inhibitory': 0.1,
              'Sigma inhibitory': 0.05,
              'Number excitational cells': 8,
              'Number inhibitory cells': 2,
              'Weight Stimulus': -300.,
              'Radius stimuls': 0.2,
              'Sigma Stimuls': 0.1,
              'Background rate': 45000.,
              'Time before stimulation': 500.,
              'Time of stimulation': 200.,
              'Time after Stimulation': 500.,
              }

pd.set_option('display.max_columns', 500)

nest.SetKernelStatus({"resolution": 0.1, "print_time": True, "overwrite_files": True})

nest.CopyModel('iaf_psc_alpha', 'exci')
nest.CopyModel('iaf_psc_alpha', 'inhi')
nest.CopyModel('static_synapse', 'exc', {'weight': 5.0})
nest.CopyModel('static_synapse', 'inh', {'weight': -20.0})
nest.CopyModel('static_synapse', 'inh_strong', {'weight': 30.0})

l = tp.CreateLayer({'rows': cols_rows,
                    'columns': cols_rows,
                    'elements': ['exci', 8, 'inhi', 2]})

cdict_e2i = {'connection_type': 'divergent',
             'mask': {'circular': {'radius': 0.05}},
             'kernel': {'gaussian': {'p_center': 0.8, 'sigma': .025}},
             'delays': {'linear': {'c': 2.0, 'a': 0.02}},
             'sources': {'model': 'exci'},
             'targets': {'model': 'inhi'},
             'synapse_model': 'exc'}
cdict_e2e = {'connection_type': 'divergent',
             'mask': {'circular': {'radius': 0.05}},
             'kernel': {'gaussian': {'p_center': 0.8, 'sigma': .025}},
             'delays': {'linear': {'c': 2.0, 'a': 0.02}},
             'sources': {'model': 'exci'},
             'targets': {'model': 'exci'},
             'synapse_model': 'exc'}
cdict_i2e = {'connection_type': 'divergent',
             'mask': {'circular': {'radius': 0.1}},
             'kernel': {'gaussian': {'p_center': 0.8, 'sigma': 0.05}},
             'delays': {'linear': {'c': 4.0, 'a': 0.04}},
             'sources': {'model': 'inhi'},
             'targets': {'model': 'exci'},
             'synapse_model': 'inh'}
cdict_i2i = {'connection_type': 'divergent',
             'mask': {'circular': {'radius': 0.1}},
             'kernel': {'gaussian': {'p_center': 0.8, 'sigma': 0.05}},
             'delays': {'linear': {'c': 4.0, 'a': 0.04}},
             'sources': {'model': 'inhi'},
             'targets': {'model': 'inhi'},
             'synapse_model': 'inh'}

tp.ConnectLayers(l, l, cdict_e2i)
tp.ConnectLayers(l, l, cdict_e2e)
tp.ConnectLayers(l, l, cdict_i2i)
tp.ConnectLayers(l, l, cdict_i2e)


stim = tp.CreateLayer({'rows': 1,
                       'columns': 1,
                       'elements': 'poisson_generator'})
stim_i = nest.GetLeaves(stim, local_only=True)[0]
nest.SetStatus(stim_i, {'rate': 45000.0})
cdict_stim = {'connection_type': 'divergent',
              'mask': {'circular': {'radius': 2.}},
              'synapse_model': 'exc'}

tp.ConnectLayers(stim, l, cdict_stim)

stim2 = tp.CreateLayer({'rows': 1,
                       'columns': 1,
                       'elements': 'poisson_generator'})
stim2_i = nest.GetLeaves(stim2, local_only=True)[0]
nest.SetStatus(stim2_i, {'rate': 0.0})
cdict_stim2 = {'connection_type': 'divergent',
               'kernel': {'gaussian': {'p_center': 1., 'sigma': 0.1}},
               'mask': {'circular': {'radius': 0.2},
                        'anchor': [0., 0.]},
               'targets': {'model': 'exci'},
               'synapse_model': 'inh_strong'}

tp.ConnectLayers(stim2, l, cdict_stim2)


rec = nest.Create("spike_detector")
nrns = nest.GetLeaves(l, local_only=True)[0]

nest.Connect(nrns, rec)

rec_ex = tp.CreateLayer({'rows': 1,
                         'columns': 1,
                         'elements': 'spike_detector'})
cdict_rec_ex = {'connection_type': 'convergent',
                'sources': {'model': "exci"}}
tp.ConnectLayers(l, rec_ex, cdict_rec_ex)

rec_in = tp.CreateLayer({'rows': 1,
                         'columns': 1,
                         'elements': 'spike_detector'})
cdict_rec_in = {'connection_type': 'convergent',
                'sources': {'model': 'inhi'}}
tp.ConnectLayers(l, rec_in, cdict_rec_in)

nest.Simulate(500.0)
#nest.SetStatus(nrns, {'I_e': 500.0})
nest.SetStatus(stim2_i, {'rate': 600000.0})
nest.Simulate(200.0)
#nest.SetStatus(nrns, {'I_e': 0.0})
nest.SetStatus(stim2_i, {'rate': 0.0})
nest.Simulate(200.0)

rec_ex_true = nest.GetLeaves(rec_ex, local_only=True)[0]
rec_in_true = nest.GetLeaves(rec_in, local_only=True)[0]

events_ex = nest.GetStatus(rec_ex_true, "events")[0]
events_in = nest.GetStatus(rec_in_true, "events")[0]
df_ex = makePandas(events_ex, tp.FindCenterElement(l)[0])
df_in = makePandas(events_in, tp.FindCenterElement(l)[0])

#print('Calculating Fano Factor (may take longer)')
#fano, tFano = fanoFactor(df_ex, gridSize=gridSize, tMin=0., tMax=550., tStep=25.)
#plt.plot(tFano, fano)
#plt.title('Fano over time')
#plt.xlabel('Time (ms)')
#plt.ylabel('Fano Factor')
#plt.show()
#print('Finished calculating Fano Factor.')


df_ex_before_stim = df_ex[df_ex['Time']<=500.]
visualization(df_ex_before_stim, 'Excitatory Neurons before stimulation, Excitation of Excitatory Neurons')
plt.show()
df_ex_while_stim = df_ex[df_ex['Time']>500.]
df_ex_while_stim = df_ex_while_stim[df_ex_while_stim['Time']<=700.]
visualization(df_ex_while_stim, 'Excitatory Neurons while stimulation, Excitation of Excitatory Neurons')
plt.show()
df_ex_after_stim = df_ex[df_ex['Time']>700.]
visualization(df_ex_after_stim, 'Excitatory Neurons after stimulation, Excitation of Excitatory Neurons')
plt.show()
df_in_before_stim = df_in[df_in['Time']<=500.]
visualization(df_in_before_stim, 'Inhibitory Neurons before stimulation, Excitation of Excitatory Neurons')
plt.show()
df_in_while_stim = df_in[df_in['Time']>500.]
df_in_while_stim = df_in_while_stim[df_in_while_stim['Time']<=700.]
visualization(df_ex_while_stim, 'Inhibitory Neurons while stimulation, Excitation of Excitatory Neurons')
plt.show()
df_in_after_stim = df_in[df_in['Time']>700.]
visualization(df_in_after_stim, 'Inhibitory Neurons after stimulation, Excitation of Excitatory Neurons')
plt.show()

wvisualization(df_ex, 'Excitatory Neurons, Excitation of Excitatory Neurons')
plt.savefig(all_folder+'visu_exci_exci.pdf', dpi=300)
visualization(df_in, 'Inhibitory Neurons, Excitation of Excitatory Neurons')
plt.savefig(all_folder+'visu_inhi_exci.pdf', dpi=300)

distanceFiringRate(events_ex, tp.FindCenterElement(l)[0], neurons_per_gridpoint=8, title="Excitatory", gridSize=gridSize)
plt.savefig(all_folder+'dist_exci_exci.pdf', dpi=300)
distanceFiringRate(events_in, tp.FindCenterElement(l)[0], neurons_per_gridpoint=2, title='Inhibitory', gridSize=gridSize)
plt.savefig(all_folder+'dist_inhi_exci.pdf', dpi=300)


a,b = distance(events_ex, 0, 0.1, tp.FindCenterElement(l)[0])
a1,b1 = distance(events_ex, 0.1, 0.2, tp.FindCenterElement(l)[0])
c,d = distance(events_in, 0, 0.1, tp.FindCenterElement(l)[0])
c1,d1 = distance(events_in, 0.1, 0.2, tp.FindCenterElement(l)[0])
raster_plot(senders=a, timeS=b, title="Distance from 0 to 0.1 (Excitatory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_01_raster_exci_exci.pdf', dpi=300)
raster_plot(senders=c, timeS=d, title="Distance from 0 to 0.1 (Inhbitory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_01_raster_inhi_exci.pdf', dpi=300)
raster_plot(senders=a1, timeS=b1, title="Distance from 0.1 to 0.2 (Excitatory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_02_raster_exci_exci.pdf', dpi=300)
raster_plot(senders=c1, timeS=d1, title="Distance from 0.1 to 0.2 (Inhbitory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_02_raster_inhi_exci.pdf', dpi=300)
a2,b2 = distance(events_ex, 0.2, 0.3, tp.FindCenterElement(l)[0])
a3,b3 = distance(events_ex, 0.3, 0.4, tp.FindCenterElement(l)[0])
c2,d2 = distance(events_in, 0.2, 0.3, tp.FindCenterElement(l)[0])
c3,d3 = distance(events_in, 0.3, 0.4, tp.FindCenterElement(l)[0])
raster_plot(senders=a2, timeS=b2, title="Distance from 0.2 to 0.3 (Excitatory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_03_raster_exci_exci.pdf', dpi=300)
raster_plot(senders=a3, timeS=b3, title="Distance from 0.3 to 0.4 (Excitatory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_04_raster_exci_exci.pdf', dpi=300)
raster_plot(senders=c2, timeS=d2, title="Distance from 0.2 to 0.3 (Inhibitory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_03_raster_inhi_exci.pdf', dpi=300)
raster_plot(senders=c3, timeS=d3, title="Distance from 0.3 to 0.4 (Inhibitory Neurons/ Excitation of Excitatory Neurons)", gridSize=gridSize)
plt.savefig(all_folder+'dist_04_raster_inhi_exci.pdf', dpi=300)

#print(rec)
#nest.raster_plot.from_device(rec_ex_true, hist=True)
#nest.raster_plot.show()
#nest.raster_plot.from_device(rec_in_true, hist=True)
#nest.raster_plot.show()

#fig = tp.PlotLayer(l)
#ctr = tp.FindCenterElement(l)
#tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
#plt.show()
#fig = tp.PlotLayer(l)
#ctr = tp.FindCenterElement(l)
#tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')
raster_plot(eventSenders=events_ex, title='Excitatory Neurons', gridSize=gridSize)
plt.savefig(all_folder+'raster_exci_exci.pdf', dpi=300)
raster_plot(eventSenders=events_in, title='Inhibitory Neurons', gridSize=gridSize)
plt.savefig(all_folder+'raster_inhi_exci.pdf', dpi=300)
#print(rec)
#nest.raster_plot.from_device(rec_ex_true, hist=True)
#nest.raster_plot.show()
#nest.raster_plot.from_device(rec_in_true, hist=True)
#nest.raster_plot.show()




#fig = tp.PlotLayer(l)
#ctr = tp.FindCenterElement(l)
#tp.PlotTargets(ctr, l, fig=fig, mask=cdict_e2e['mask'], kernel=cdict_e2e['kernel'], src_size=250, tgt_color='red', tgt_size=20, kernel_color='green', syn_type='exc')
#plt.show()
#fig = tp.PlotLayer(l)
#ctr = tp.FindCenterElement(l)
#tp.PlotTargets(ctr, l, fig=fig, mask=cdict_i2e['mask'], kernel=cdict_i2e['kernel'], src_size=200, tgt_color='red', tgt_size=10, kernel_color='green')
#plt.show()

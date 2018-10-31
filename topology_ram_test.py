import nest
import gc
import nest.topology as tp
import numpy as np


# Short Class for creating network, to see if ram gets deallocated
class RandomBalancedNetwork:
    def __init__(self):
        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": 0.1, "print_time": True,
                              "overwrite_files": True,
                              "local_num_threads": 8
                             })
        nest.CopyModel('iaf_psc_alpha', 'exci')
        nest.CopyModel('iaf_psc_alpha', 'inhi')
        nest.CopyModel('static_synapse', 'exc', {'weight': 5.0})
        nest.CopyModel('static_synapse', 'inh', {'weight': -5.0})
        self.l = tp.CreateLayer({'rows': 90,
                                 'columns': 90,
                                 'elements': ['exci', 5,
                                              'inhi', 5],
                                 'edge_wrap': False})
        cdict = {'connection_type': 'divergent',
                     'mask': {'circular': {'radius': 0.2}},
                     'kernel': {'gaussian': {'p_center': 0.8, 'sigma': 0.075}},
                     'delays': {'linear': {'c': 2.0, 'a': 0.02}},
                     'sources': {'model': 'exci'},
                     'targets': {'model': 'inhi'},
                     'synapse_model': 'exc'}
        tp.ConnectLayers(self.l, self.l, cdict)
        self.rec_ex = tp.CreateLayer({'rows': 1,
                                      'columns': 1,
                                      'elements': 'spike_detector'})
        cdict_rec_ex = {'connection_type': 'convergent',
                        'sources': {'model': "exci"}}
        tp.ConnectLayers(self.l, self.rec_ex, cdict_rec_ex)
        # Background stimulation
        stim = tp.CreateLayer({'rows': 1,
                               'columns': 1,
                               'elements': 'poisson_generator'})
        stim_i = nest.GetLeaves(stim, local_only=True)[0]
        nest.SetStatus(stim_i, {'rate': 30000.})
        background_stim_dict = {'connection_type': 'divergent',
                                'mask': {'grid': {'rows': 90,
                                                  'columns': 90}},
                                'synapse_model': 'exc'}
        tp.ConnectLayers(stim, self.l, background_stim_dict)
        nest.Simulate(2000.)
        rec_ex_true = nest.GetLeaves(self.rec_ex, local_only=True)[0]
        self.events_ex = nest.GetStatus(rec_ex_true, "events")[0]

for i in range(100):
    print("Run: " + str(i))
    sim = RandomBalancedNetwork()

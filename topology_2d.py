import pandas as pd
import os
import topology_2d_helper as magic
import sys
import general_helper as gen
import matplotlib
matplotlib.use('Agg')


base_folder = 'figures/'
date_folder = 'tsodyks_analysis/'
sub_folder = '18-12-07/'
all_folder = base_folder + date_folder + sub_folder
pd.set_option('display.max_columns', 500)

parameters = {'Name': 'tsodyks',
              'Columns': 40,
              'Rows': 40,
              'Radius excitational': 0.2,
              'Sigma excitational': 0.1,
              'Radius inhibitory': 0.2,
              'Sigma inhibitory': 0.1,
              'Jee': 3.,
              'Jii': 1.,
              'Jei': 3.,
              'Jie': -3.,
              'Jee Connectivity': 0.4,
              'Jii Connectivity': 0.4,
              'Jei Connectivity': 0.4,
              'Jie Connectivity': 0.4,
              'Number excitational cells': 5,
              'Number inhibitory cells': 5,
              'Weight Stimulus': -10., # old: -3000.
              'Radius stimulus': 0.1,
              'Sigma Stimulus': 0.05,
              'e2e delay': 1.0,
              'e2i delay': 1.0,
              'i2e delay': 1.0,
              'i2i delay': 1.0,
              'delay growth multiplier': 2,
              'Stimulus rate': 40000.,
              'Background rate excitatory': 25000.,
              'Background rate inhibitory': 15000.,
              'Time before stimulation': 500.,
              'Time of stimulation': 100.,
              'Time after Stimulation': 500.,
              }

"""
Starter
"""
def tsodyksAnalysis():
    curr_folder = all_folder + "40x40/test_higher_synpatic_delay"
    magic.makeDir(curr_folder)
    magic.tsodyks_analysis(parameters, curr_folder)

def test():
    curr_folder = base_folder + sub_folder + "test1_inh_decrease_middle"
    magic.makeDir(curr_folder+'/excitatory_neurons')
    magic.makeDir(curr_folder+'/inhibitory_neurons')
    magic.simulationAndAnalysis(parameters, curr_folder)

path = '/home/adrossel/Magic/pickle/'
current_subfolder = 'jii_3_edge_wrap'
magic.simulate_and_dump(parameters, path+current_subfolder)

"""
Run Program
"""
# magic.average_firng_rates(parameters)
# magic.tsodyks_analysis_quiver(parameters)
# tsodyksAnalysis()
# test()

def automation():
    radius_inhibs = [0.1, 0.2, 0.2, 0.2]
    sigma_inhibs = [0.075, 0.1, 0.1, 0.75]
    radius_excis = [0.1, 0.2, 0.1, 0.1]
    sigma_excis = [0.05, 0.1, 0.05, 0.5]
    weight_inhis = [-5., -5.0, -4.0, -20.0]
    weight_excis = [5., 5.0, 4.0, 20.0]
    sim_folder = base_folder + date_folder
    for cols_rows in [30, 40, 45 , 50]:
        col_folder = sim_folder + '/colsRows_' + str(cols_rows)
        gen.create_folder(col_folder)
        parameters['Columns'] = cols_rows
        parameters['Rows'] = cols_rows
        for background_rate in [30000., 35000.]:
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
                magic.simulationAndAnalysis(parameters, curr_folder)


def testSimulation():
    radius_excis = [0.1, 0.2, 0.1, 0.2, 0.1]
    radius_inhibs = [0.1, 0.2, 0.2, 0.1, 0.2]
    sigma_excis = [0.05, 0.1, 0.05, 0.1, 0.5]
    sigma_inhibs = [0.075, 0.1, 0.1, 0.05, 0.75]
    weight_inhis = [-5., -5.0, -4.0, -10.0, -20.0]
    weight_excis = [5., 5.0, 4.0, 10.0, 20.0]
    sim_folder = base_folder + date_folder
    for cols_rows in [80]:
        col_folder = sim_folder + '/colsRows_' + str(cols_rows)
        gen.create_folder(col_folder)
        parameters['Columns'] = cols_rows
        parameters['Rows'] = cols_rows
        for background_rate in [25000., 30000., 35000.]:
            background_folder = col_folder + '/background_' + str(background_rate)
            gen.create_folder(background_folder)
            parameters['Background rate'] = background_rate
            for radius_inhib, sigma_inhib, radius_exci, sigma_exci, weight_inhi, weight_exci in zip(radius_inhibs,
                                                sigma_inhibs, radius_excis, sigma_excis, weight_inhis, weight_excis):
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
                print('Loading parameters in Simulation.')
                magic.randBalNetwork(parameters)


def bash():
    print(str(sys.argv))
    parameters['Columns'] = int(sys.argv[1])
    parameters['Rows'] = int(sys.argv[2])
    parameters['Background Rate'] = float(sys.argv[3])
    parameters['Radius excitational'] = float(sys.argv[4])
    parameters['Radius inhibitory'] = float(sys.argv[5])
    parameters['Sigma excitational'] = float(sys.argv[6])
    parameters['Sigma inhibitory'] = float(sys.argv[7])
    parameters['Inhibitory Weight'] = float(sys.argv[8])
    parameters['Excitational Weight'] = float(sys.argv[9])
    parameters['Weight Stimulus'] = float(sys.argv[10])
    sim_folder = base_folder + date_folder
    col_folder = sim_folder + '/colsRows_' + str(parameters['Columns'])
    gen.create_folder(col_folder)
    background_folder = col_folder + '/background_' + str(parameters['Background Rate'])
    gen.create_folder(background_folder)
    name = "inh_(" + str(parameters['Radius inhibitory']) + "," + str(parameters['Sigma inhibitory']) + ")_exci_(" + str(parameters['Radius excitational']) + "," + str(parameters['Sigma excitational']) + "_weight_inhi_" + str(parameters['Inhibitory Weight']) + "weight_exci_" + str(parameters['Excitational Weight']) + ")"
    curr_folder = background_folder + '/' + name
    parameters['Name'] = name
    if not os.path.exists(curr_folder):
        os.mkdir(curr_folder)
        os.mkdir(curr_folder+'/excitatory_neurons')
        os.mkdir(curr_folder+'/inhibitory_neurons')
    magic.simulationAndAnalysis(parameters, curr_folder)




#automation()
#testSimulation()
#bash()

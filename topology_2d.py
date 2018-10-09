import sys
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
date_folder = '18_10_09/'
sub_folder = 'tests/'
all_folder = base_folder + date_folder + sub_folder
pd.set_option('display.max_columns', 500)

parameters = {'Name': 'inhibition',
              'Columns': 80,
              'Rows': 80,
              'Excitational Weight': 5.0,
              'Radius excitational': 0.05,
              'Sigma excitational': 0.025,
              'Inhibitory Weight': -4.0,
              'Radius inhibitory': 0.075,
              'Sigma inhibitory': 0.0375,
              'Number excitational cells': 5,
              'Number inhibitory cells': 5,
              'Weight Stimulus': -3000.,
              'Radius stimulus': 0.1,
              'Sigma Stimulus': 0.05,
              'Stimulus rate': 40000.,
              'Background rate': 19000.,
              'Time before stimulation': 1000.,
              'Time of stimulation': 500.,
              'Time after Stimulation': 1000.,
              }

def test():
    curr_folder = all_folder + "new"
    if not os.path.exists(curr_folder):
        os.mkdir(curr_folder)
        os.mkdir(curr_folder+'/excitatory_neurons')
        os.mkdir(curr_folder+'/inhibitory_neurons')
    magic.simulationAndAnalysis(parameters, curr_folder)


test()

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
    radius_excis =  [0.1, 0.2, 0.1, 0.2, 0.1]
    radius_inhibs = [0.1, 0.2, 0.2, 0.1, 0.2]
    sigma_excis =   [0.05, 0.1, 0.05, 0.1, 0.5]
    sigma_inhibs =  [0.075, 0.1, 0.1, 0.05, 0.75]
    weight_inhis =  [-5., -5.0, -4.0, -10.0, -20.0]
    weight_excis =  [5., 5.0, 4.0, 10.0, 20.0]
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

import json
import numpy as numpy
import matplotlib.pyplot as plt


placebo_RR50_array_file_name = 'placebo_RR50_array'
placebo_MPC_array_file_name = 'placebo_MPC_array'
drug_RR50_array_file_name = 'drug_RR50_array'
drug_MPC_array_file_name = 'drug_MPC_array'

def load_endpoint_from_trial(endpoint_file_name, trial_ID_num):

    with open( os.getcwd() + '/' + endpoint_file_name + '.json', 'r')


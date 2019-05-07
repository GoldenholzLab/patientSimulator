import json
import numpy as np
import matplotlib.pyplot as plt
import os
import textwrap
import sys


def load_data_point_from_trial(endpoint_file_name, trial_ID_num):
    '''

    Inputs:

        1) endpoint_file_name:

            (string) - the file name of the 

    Outputs:

    '''

    with open( os.getcwd() + '/' + endpoint_file_name + '.json', 'r') as json_file:

        endpoint_array = np.array(json.load(json_file))
    
    return endpoint_array[trial_ID_num]


def load_data_points_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                              drug_RR50_array_file_name,    drug_MPC_array_file_name,
                              RR50_p_values_filename,       MPC_p_values_filename,
                              trial_ID_num,                 num_patients_per_arm_file_name):
    '''

    Inputs:

        1)

    Outputs:

    '''

    placebo_RR50 = load_data_point_from_trial(placebo_RR50_array_file_name, trial_ID_num)
    placebo_MPC  = load_data_point_from_trial(placebo_MPC_array_file_name,  trial_ID_num)
    drug_RR50    = load_data_point_from_trial(drug_RR50_array_file_name,    trial_ID_num)
    drug_MPC     = load_data_point_from_trial(drug_MPC_array_file_name,     trial_ID_num)
    RR50_p_value = load_data_point_from_trial(RR50_p_values_filename,       trial_ID_num)
    MPC_p_value  = load_data_point_from_trial(MPC_p_values_filename,        trial_ID_num)

    with open( os.getcwd() + '/' + num_patients_per_arm_file_name + '.json', 'r') as json_file:

        num_patients_per_arm = json.load(json_file)

    return [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, num_patients_per_arm]


def plot_trial_endpoints(placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, 
                         num_patients_per_arm, endpoint_decimal_round, p_value_decimal_round):

    # Plot the 50% Responder Rate in the usual way that other RCT papers graph it
    [fig, (ax1, ax2)] = plt.subplots(1, 2)
    rects1 = ax1.bar(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], [drug_RR50, placebo_RR50], color = 'k')
    ax1.set_xticklabels(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], fontsize = 14)
    ax1.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax1.set_ylim([0, 80])
    ax1.set_ylabel('Percentage of 50% Responders', fontsize = 14)
    longTitle = '50% Responder Rate Results of One Trial (p' + str( np.round( RR50_p_value, p_value_decimal_round ) ) + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 40))
    ax1.set_title(formatted_title, fontsize = 14)
    for rect in rects1:
        height = rect.get_height()
        ax1.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, endpoint_decimal_round), ha = 'center', va = 'bottom', fontsize = 14)
        
    # Plot the Median Percent Change in the usual way that other RCT papers graph it
    rects2 = ax2.bar(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], 
                      [drug_MPC, placebo_MPC], color = 'k')
    ax2.set_xticklabels(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], fontsize = 14)
    ax2.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax2.set_ylim([0, 80])
    ax2.set_ylabel('Median Percentage Change', fontsize = 14)
    longTitle = 'Median Percentage Change Results of One Trial (p' + str( np.round( MPC_p_value, p_value_decimal_round ) ) + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 45))
    ax2.set_title(formatted_title, fontsize = 14)
    for rect in rects2:
        height = rect.get_height()
        ax2.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, endpoint_decimal_round), ha = 'center', va = 'bottom', fontsize = 14)
    
    left  = -0.4  # the left side of the subplots of the figure
    right = 1.3    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.2   # the amount of height reserved for white space between subplots
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    fig.savefig(fname = os.getcwd() + '/Romero-fig5.png', dpi = 600, bbox_inches = 'tight')


if(__name__ == '__main__'):

    trial_ID_num = int(sys.argv[1])

    placebo_RR50_array_file_name = 'placebo_RR50_array'
    placebo_MPC_array_file_name = 'placebo_MPC_array'
    drug_RR50_array_file_name = 'drug_RR50_array'
    drug_MPC_array_file_name = 'drug_MPC_array'
    RR50_p_values_filename = 'RR50_p_values'
    MPC_p_values_filename = 'MPC_p_values'
    num_patients_per_arm_file_name = 'num_patients_per_arm'

    endpoint_decimal_round = 1
    p_value_decimal_round = 2

    [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, num_patients_per_arm] = \
        load_data_points_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                                    drug_RR50_array_file_name,    drug_MPC_array_file_name,
                                    RR50_p_values_filename,       MPC_p_values_filename,
                                    trial_ID_num,                 num_patients_per_arm_file_name)

    plot_trial_endpoints(placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, 
                         num_patients_per_arm, endpoint_decimal_round, p_value_decimal_round)
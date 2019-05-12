import json
import numpy as np
import matplotlib.pyplot as plt
import os
import textwrap
import sys

'''

This script plots figure 5 from the paper associated with this repository. It assumes the existence of all the 

needed JSON files containing the relevant information in the same folder as this Python script. If those don;t

exist, generate them with the data_generation.py script. It is recommended to run this from the command line

with the following command:


                    $ python plot_trial.py <trial_ID_number>

The trial_ID_number parameter is the ID number of the patient to plotted with respect to the entire synthetic

patient population used to generate figures 1 and 2 of the paper associated with this repository. Assuming nothing

was changed in the data_generation.py script, the total number of synthetic patients should be exactly 10,000 patients.

'''


def load_data_point_from_trial(data_point_file_name, trial_ID_num):
    '''

    This fucntion loads the data point for one data point that corresponds to a trial

    from a JSON file containing an array of data points, with each array element being

    the same data point from a different trial.

    Inputs:

        1) data_point_file_name:

            (string) - the file name of the json file that stores the array which has

                       the relevant data point of the trial to be plotted
        
        2) trial_ID_num:

            (int) - the number of the trial to be plotted

    Outputs:

    '''

    # locate the JSON file containing the relevant array and open it up for reading
    with open( os.getcwd() + '/' + data_point_file_name + '.json', 'r') as json_file:

        # load the entire array
        data_point_array = np.array(json.load(json_file))
    
    # return only the data point that corresponds to the trial specified by the trial ID number
    return data_point_array[trial_ID_num]


def load_data_points_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                                drug_RR50_array_file_name,    drug_MPC_array_file_name,
                                RR50_p_values_filename,       MPC_p_values_filename,
                                trial_ID_num,                 num_patients_per_arm_file_name):
    '''

    Inputs:

        1) placebo_RR50_array_file_name:

            (string) - the file name of the JSON file containing the array of 50% responder rates

                       from the placebo arm of every trial
        
        2) placebo_MPC_array_file_name:

            (string) - the file name of the JSON file containing the array of median percent changes

                       from the placebo arm of every trial
        
        3) drug_RR50_array_file_name:

            (string) - the file name of the JSON file containing the array of 50% responder rates

                       from the drug arm of every trial
        
        4) drug_MPC_array_file_name:

            (string) - the file name of the JSON file containing the array of median percent changes

                       from the drug arm of every trial
        
        5) RR50_p_values_filename:

            (string) - the file name of the JSON file containing the array of p-values for the 50%
                       
                       responder rate from all trials
        
        6) MPC_p_values_filename:

            (string) - the file name of the JSON file containing the array of p-values for the
                       
                       median percent change from all trials
        
        7) trial_ID_num:

            (int) - the ID number of the trial to be plotted
        
        8) num_patients_per_arm_file_name:

            (string) - the file name of the JSON file containing the number of patients for both the 
                       
                       placebo arm and the drug arm

    Outputs:

        1) placebo_RR50:

            (float) - the 50% responder rate for the placebo arm of the trial to be plotted

        2) placebo_MPC:

            (float) - the median percent change for the placebo arm of the trial to be plotted

        3) drug_RR50:

            (float) - the 50% responder rate for the drug arm of the trial to be plotted

        4) drug_MPC:

            (float) - the median percent change for the drug arm of the trial to be plotted

        5) RR50_p_value:

            (float) - the p-value of the 50% responder rate of the trial to be plotted (Fisher Exact test)

        6) MPC_p_value:

            (float) - the p-value of the median percent change of the trial to be plotted (Wilcoxon Signed Rank Test)

        7) num_patients_per_arm:
            
            (int) - the number of patients generated for both arms in all trials

    '''

    # load all the endpoint responses from their respective JSON files
    placebo_RR50 = load_data_point_from_trial(placebo_RR50_array_file_name, trial_ID_num)
    placebo_MPC  = load_data_point_from_trial(placebo_MPC_array_file_name,  trial_ID_num)
    drug_RR50    = load_data_point_from_trial(drug_RR50_array_file_name,    trial_ID_num)
    drug_MPC     = load_data_point_from_trial(drug_MPC_array_file_name,     trial_ID_num)

    # load the p-values from their respective JSON files
    RR50_p_value = load_data_point_from_trial(RR50_p_values_filename,       trial_ID_num)
    MPC_p_value  = load_data_point_from_trial(MPC_p_values_filename,        trial_ID_num)

    # manually load the number of patients per arm from a JSON file
    with open( os.getcwd() + '/' + num_patients_per_arm_file_name + '.json', 'r') as json_file:

        num_patients_per_arm = json.load(json_file)

    return [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, num_patients_per_arm]


def plot_trial_endpoints(placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, 
                         num_patients_per_arm, endpoint_decimal_round, p_value_decimal_round):
    '''

    Inputs:

        1) placebo_RR50:

            (float) - the 50% responder rate for the placebo arm of the trial to be plotted

        2) placebo_MPC:

            (float) - the median percent change for the placebo arm of the trial to be plotted

        3) drug_RR50:

            (float) - the 50% responder rate for the drug arm of the trial to be plotted

        4) drug_MPC:

            (float) - the median percent change for the drug arm of the trial to be plotted

        5) RR50_p_value:

            (float) - the p-value of the 50% responder rate of the trial to be plotted (Fisher Exact test)

        6) MPC_p_value:

            (float) - the p-value of the median percent change of the trial to be plotted (Wilcoxon Signed Rank Test)
        
        7) num_patients_per_arm:
            
            (int) - the number of patients generated for both arms in all trials
        
        8) endpoint_decimal_round:

            (int) - the number of decimal places to round the endpoint responses too
        
        9) p_value_decimal_round:

            (int) - the number of decimal places to round the p-values too

    Outputs:

        Technically None

    '''

    # Plot the 50% Responder Rate in the usual way that other RCT papers graph it
    [fig, (ax1, ax2)] = plt.subplots(1, 2)
    rects1 = ax1.bar(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], [drug_RR50, placebo_RR50], color = 'k')
    ax1.set_xticklabels(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], fontsize = 14)
    ax1.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax1.set_ylim([0, 80])
    ax1.set_ylabel('Percentage of 50% Responders', fontsize = 14)
    longTitle = '50% Responder Rate Results of One Trial (p = ' + str( np.round( RR50_p_value, p_value_decimal_round ) ) + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 40))
    ax1.set_title(formatted_title, fontsize = 14)
    for rect in rects1:
        height = rect.get_height()
        ax1.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%r'%np.round(height, endpoint_decimal_round), ha = 'center', va = 'bottom', fontsize = 14)
    
    # Plot the Median Percent Change in the usual way that other RCT papers graph it
    rects2 = ax2.bar(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], 
                      [drug_MPC, placebo_MPC], color = 'k')
    ax2.set_xticklabels(['drug (' + str( num_patients_per_arm ) + ' patients)', 'placebo (' + str( num_patients_per_arm ) + ' patients)'], fontsize = 14)
    ax2.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax2.set_ylim([0, 80])
    ax2.set_ylabel('Median Percentage Change', fontsize = 14)
    longTitle = 'Median Percentage Change Results of One Trial (p = ' + str( np.round( MPC_p_value, p_value_decimal_round ) ) + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 45))
    ax2.set_title(formatted_title, fontsize = 14)
    for rect in rects2:
        height = rect.get_height()
        ax2.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%r'%np.round(height, endpoint_decimal_round), ha = 'center', va = 'bottom', fontsize = 14)
    
    left  = -0.4  # the left side of the subplots of the figure
    right = 1.3    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.2   # the amount of height reserved for white space between subplots
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    fig.savefig(fname = os.getcwd() + '/Romero-fig5.png', dpi = 600, bbox_inches = 'tight')


if(__name__ == '__main__'):

    # take the ID number of the trial to be plotted from the command line
    '''
    
    For current data, this should be 26

    '''
    trial_ID_num = int(sys.argv[1])

    # the file names of the JSON files that all the data points of the JSON files are stored in
    placebo_RR50_array_file_name = 'placebo_RR50_array'
    placebo_MPC_array_file_name = 'placebo_MPC_array'
    drug_RR50_array_file_name = 'drug_RR50_array'
    drug_MPC_array_file_name = 'drug_MPC_array'
    RR50_p_values_filename = 'RR50_p_values'
    MPC_p_values_filename = 'MPC_p_values'
    num_patients_per_arm_file_name = 'num_patients_per_arm'

    # the number of decimal places to round the endpoint responses too
    endpoint_decimal_round = 1

    # the number of decimal places to round the p-values too
    p_value_decimal_round = 3

    # load all the relevant data points
    [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, num_patients_per_arm] = \
        load_data_points_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                                    drug_RR50_array_file_name,    drug_MPC_array_file_name,
                                    RR50_p_values_filename,       MPC_p_values_filename,
                                    trial_ID_num,                 num_patients_per_arm_file_name)

    # plot the tiral using the loaded data points
    plot_trial_endpoints(placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, RR50_p_value, MPC_p_value, 
                         num_patients_per_arm, endpoint_decimal_round, p_value_decimal_round)
import json
import numpy as numpy
import matplotlib.pyplot as plt


def load_endpoint_from_trial(endpoint_file_name, trial_ID_num):

    with open( os.getcwd() + '/' + endpoint_file_name + '.json', 'r') as json_file:

        endpoint_array = np.array(json.load(json_file))
    
    return endpoint_array[trial_ID_num]


def load_endpoints_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                              drug_RR50_array_file_name,    drug_MPC_array_file_name, 
                              trial_ID_num,                 num_patients_per_arm_file_name):

    placebo_RR50 = load_endpoint_from_trial(placebo_RR50_array_file_name, trial_ID_num)
    placebo_MPC  = load_endpoint_from_trial(placebo_MPC_array_file_name , trial_ID_num)
    drug_RR50    = load_endpoint_from_trial(drug_RR50_array_file_name   , trial_ID_num)
    drug_MPC     = load_endpoint_from_trial(drug_MPC_array_file_name    , trial_ID_num)

    with open( os.getcwd(() + '/' + num_patients_per_arm_file_name + '.json', 'r') as json_file:

        num_patients_per_arm = json.load(json_file)

    return [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, num_patients_per_arm]


def plot_trial_endpoints(placebo_RR50, placebo_MPC, drug_RR50, drug_MPC):

    [placebo_RR50, placebo_MPC, drug_RR50, drug_MPC, num_patients_per_arm] = \
        load_endpoints_from_trial(placebo_RR50_array_file_name, placebo_MPC_array_file_name, 
                                  drug_RR50_array_file_name,    drug_MPC_array_file_name, 
                                  trial_ID_num,                 num_patients_per_arm_file_name)

    

    # Plot the 50% Responder Rate in the usual way that other RCT papers graph it
    [fig, (ax1, ax2)] = plt.subplots(1, 2)
    rects1 = ax1.bar(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], [100*drugRR50Output, 100*placeboRR50Output], color = 'k')
    ax1.set_xticklabels(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], fontsize = 14)
    ax1.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax1.set_ylim([0, 80])
    ax1.set_ylabel('Percentage of 50% Responders', fontsize = 14)
    longTitle = '50% Responder Rate Results of One Trial (p' + RR50PValueString + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 40))
    ax1.set_title(formatted_title, fontsize = 14)
    for rect in rects1:
        height = rect.get_height()
        ax1.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, genValueRounding), ha = 'center', va = 'bottom', fontsize = 14)
        
    # Plot the Median Percent Change in the usual way that other RCT papers graph it
    rects2 = ax2.bar(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], 
                      [100*drugMPCOutput, 100*placeboMPCOutput], color = 'k')
    ax2.set_xticklabels(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], fontsize = 14)
    ax2.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax2.set_ylim([0, 80])
    ax2.set_ylabel('Median Percentage Change', fontsize = 14)
    longTitle = 'Median Percentage Change Results of One Trial (p' + MPCPValueString + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 45))
    ax2.set_title(formatted_title, fontsize = 14)
    for rect in rects2:
        height = rect.get_height()
        ax2.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, genValueRounding), ha = 'center', va = 'bottom', fontsize = 14)
    
    left  = -0.4  # the left side of the subplots of the figure
    right = 1.3    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.2   # the amount of height reserved for white space between subplots
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    fig.savefig(fname = image_dir + '/Romero-fig5', dpi = 600, bbox_inches = box)


placebo_RR50_array_file_name = 'placebo_RR50_array'
placebo_MPC_array_file_name = 'placebo_MPC_array'
drug_RR50_array_file_name = 'drug_RR50_array'
drug_MPC_array_file_name = 'drug_MPC_array'
num_patients_per_arm_file_name = 'num_patients_per_arm'
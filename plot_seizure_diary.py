import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import textwrap
import sys
import os
from PIL import Image
import io

'''

This script plots figure 3 from the paper associated with this repository. It assumes the existence of all the 

needed JSON files containing the relevant information in the same folder as this Python script. If those don;t

exist, generate them with the data_generation.py script. It is recommended to run this from the command line

with the following command:


                    $ python plot_patient.py <patient_ID_number>

The patient_ID_number parameter is the ID number of the patient to plotted with respect to the entire synthetic

patient population used to generate figures 1 and 2 of the paper associated with this repository. Assuming nothing

was changed in the data_generation.py script, the total number of synthetic patients should be exactly 10,000 patients.

'''

def load_patient(patient_population_weekly_seizure_counts_file_name, patient_num_id):
    '''

    This function extracts one patient from the JSON file containing the synthetic patient population, and returns

    their weekly seizure diary for plotting.

    Inputs:

        1) patient_population_weekly_seizure_counts_file_name:

            (string) - the file name of the JSON which stores the weekly seizure diaries of the synthetic patient

                       population
        
        2) patient_num_id:

            (int) - the number of the patient to be plotted from the synthetic patient patient population

    Outputs:

        1) patient_weekly_seizure_counts:

            (1D Numpy array) - the weekly seizure diary of the patient to be plotted

    '''

    # locate the JSON file containing the synthetic patient population
    with open(os.getcwd() + '/' + patient_population_weekly_seizure_counts_file_name + '.json', 'r') as json_file:

        # read the entire population into a List of seizure diaries
        patient_population_weekly_seizure_counts = json.load(json_file)

    # extract one seizure diary whose ID number is specififed by the user of this Python script
    patient_weekly_seizure_counts = patient_population_weekly_seizure_counts[patient_num_id]

    return patient_weekly_seizure_counts


def plot_patient(patient_weekly_seizure_counts):
    '''

    This function takes one weekly seizure diary and plots it

    Inputs:

        1) patient_weekly_seizure_counts:

            (1D Numpy array) - the weekly seizure diary of the patient to be plotted

    Outputs:

        Technically None

    '''

    fig = plt.figure()
    ax = fig.gca()
    numWeeks = len( patient_weekly_seizure_counts )
    plt.bar( np.arange( numWeeks ), patient_weekly_seizure_counts, color='k')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Weeks', fontsize = 14)
    plt.ylabel('weekly seizure counts', fontsize = 14)
    long_title = 'One patient\'s weekly seizure counts over ' + str(numWeeks) + ' weeks'
    formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
    plt.title(formatted_title, fontsize = 14)

    # save as picture as TIFF file
    png1 = io.BytesIO()
    fig.savefig(png1, dpi = 600, bbox_inches = 'tight', format='png')
    png2 = Image.open(png1)
    png2.save('Romero-fig3.tiff')
    png1.close()


if(__name__ == '__main__'):
    
    # the file name of the JSON file containing the entire synthetic patient population
    '''
    
    For current data, this should be 240

    '''
    patient_population_weekly_seizure_counts_file_name = 'synthetic_patient_population'

    # take the ID number of the patient to be plotted from the command line
    patient_num_id = int(sys.argv[1])

    # extract and plot the weekly seizure diary of the synthetic patient with that specific ID number
    patient_weekly_seizure_counts = load_patient(patient_population_weekly_seizure_counts_file_name, patient_num_id)
    plot_patient(patient_weekly_seizure_counts)


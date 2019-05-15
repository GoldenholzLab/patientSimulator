# patientSimulator
synthetic patient seizure diary simulator and RCT simulator

WHAT IS THIS REPOSITORY FOR?:

    This repository is meant for providing access to the code which implemented the mathematical model presented in the 
    Seizure Diary Simulator paper.

In order to generate one patient:

    The data_generation.py script has a function called generate_patient(). The documentation within the function says more,
    but essentially, you have to eenter the NV model parameters along with the lengths of the baseline and testing period 
    (the baseline period can be made to have a length of zero if need be).

In order to generate the graphs associated with the Seizure Diary Simulator paper:

    Run these scripts from the command line with the following commands in this order:

        $ python data_generation.py
        $ python plotting_data.py
        $ python plot_seizure_diary.py <patient_ID_number>
        $ python plot_trial.py <trial_ID_number>
    
    The documentation within the third and fourth scripts says more about what their command line arguments are, but those are 
    essentially the ID numbers of the patient and trial plotted in figures 3 and 5, respectively. Assuming everything executes
    correctly, the entire process should take less than 10 minutes and should generate several JSON files as well as the
    relevant figures. It is important that these scripts are all in the same folder when they are used to generate figures from
    scratch: the data_generation.py stores the simulated data within JSON files, and the other scripts read the JSON files
    in order to plot that data. Both scripts assume that the JSON files are in the same folder as the currently running script.

data_generation.py

    This script generates all the JSON files which contain the necessary information for creating the figures
    from the paper.

plotting_data.py

    This script plots figures 1, 2, and 4 from the paper, assuming the relevant JSON files are already in the same folder. This script takes
    no command line arguments.

plot_seizure_diary.py

    This script plots figure 3 from the paper, assuming the relevant JSON files are already in the same folder. It takes one command line
    argument: the ID number of the patient whose seizure diary is to be plotted. This ID number is with respect to all of the other patients
    whose data were used to create figures 1 and 2 from the paper.

plot_trial.py

    This script plots figure 5 from the paper, assuming the relevant JSON files are already in the same folder. It takes one command line
    argument: the ID number of the trial to be plotted. This ID number is with respect to all of the other trials whose data were used 
    to create figure 4 from the paper.

param_expl.py

    This script provides the code for doing one randomly initialized search for optimal NV model parameters 
    over the parameter space of the group level effects for NV model. It puts the final result of its 
    search into a small text file.

param_expl_wrap.sh

    This is a shell wrapper script for submitting the param_expl.py script as a batch job to the Harvard Medical School O2 computing cluster.

param_expl_wrap_coordinator.sh

    This is a shell script which submits many searches over the parameter space of the NV model to the Harvard Medical School O2 computing cluster.

output_searcher.py

    This script searches among the text files outputted by param_expl.py to find the most optimal set of parameters from amongst the set of 
    optimized parameters.

# patientSimulator
synthetic patient seizure diary simulator and RCT simulator

WHAT IS THIS REPOSITORY FOR?:


In order to generate one patient:


In order to generate one trial:


In order to generate the graphs associated with the [NAME HERE] paper:


testingTrialResults.py

    This is the main script which generates the graphs in the paper associated with this repository. It has 
    dependencies on TrialArmClass.py, PatientClass.py, seizureMods.py, RR50Test.py, and MPCTest.py

TrialArmClass.py

    This script provides the code for generating clinical trials of synthetic patients according to NV model. 
    It has dependencies on the PatientClass.py and seizureMods.py scripts.

PatientClass.py

    This script provides the code for generating one individual synthetic patient according to NV model. It 
    has a dependency on seizureMods.py.

seizureMods.py

    This script provides utility functions for applying the drug effect to synthetic patients as well as
    making sure patients comply with the RCT parameters if they are generated within the context of a trial.

RR50Test.py

    This script provides functions which evaluate the RR50 (50% responder rate) as well as the p-value (statistical significance)
    associated with the RR50, as according to the Fisher Exact test.

MPCTest.py

    This script provides functions which evaluate the MPC (median percent change) as well as the p-value (statistical significance)
    associated with the MPC, as according to the Wilcoxon Signed Rank Test.

param_expl.py

    This script provides the code for doing one randomly initialized search for optimal NV model parameters 
    over the parameter space of the group level effects for NV model. It puts the final result of its 
    search into a small text file.

param_expl_wrap.sh

    This is a shell wrapper script for submitting the param_expl.py script as batch job to the Harvard Medical School O2 computing cluster.

param_expl_wrap_coordinator.sh

    This is a shell script which submits many searches over the parameter space of the NV model to the Harvard Medical School O2 computing cluster.

output_searcher.py

    This script searches among the text files outputted by param_expl.py to find the most optimal set of parameters from amongst the set of 
    optimized parameters.

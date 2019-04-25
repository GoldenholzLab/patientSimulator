import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

'''

conversion of timescales for unrestricted patient population?

'''

def generate_patient_counts(num_intervals, num_days_per_interval):
    '''

    Generates one patient according to NV model

    Inputs:

        1) num_intervals:

            (int) - the number of intervals to generate for the seizure diary being outputted

        2) num_days_per_interval:

            (int) - the number of days that each interval containts (e.g., a week contains 7 days, whereas a month has 28 days, etc.)

    Outputs:

        1) patient_counts:

            (1D numpy array) - an array of seizure diary counts

    '''

    # NV model parameters
    #'''
    shape = 111.313
    scale = 296.728
    alpha = 296.339
    beta = 243.719
    '''
    shape= 24.143
    scale = 297.366
    alpha = 284.024
    beta = 369.628
    #'''

    # stochastically generate a new n and p for each patient
    n = np.random.gamma(shape, 1/scale)
    p = np.random.beta(alpha, beta)

    # convert to mean and overdispersion parameters
    mu = n*(1 - p)/p
    alpha = 1/n

    # initialize the array of patient seizure diary counts
    patient_counts = np.zeros(num_intervals)

    # for each interval
    for interval in range(num_intervals):

        # use the gamma-poisson mixture model to generate seizure counts via the negative binomial
        rate = np.random.gamma(num_days_per_interval/alpha, mu*alpha)
        count = np.random.poisson(rate)
        patient_counts[interval] = np.int_(count)

    return patient_counts


def generate_unrestricted_patient_population(num_patients, min_num_larger_intervals, max_num_larger_intervals, \
                                                num_days_per_smaller_interval, num_days_per_larger_interval):

    '''

    Generates patients with seizure diary lengths that are randomly generated according to integer number line. The diaries are
    randomly generated such that they can be converted to a larger timescale.

    Inputs:

        1) num_patients:

            (int) - the number of patients to generate

        2) min_num_larger_intervals:

            (int) - the minimum number of intervals to randomly generate on the larger time scale
        
        3) max_num_larger_intervals:

            (int) - the maximum number of intervals to randomly generate on the larger time scale
        
        4) num_days_per_smaller_interval:

            (int) - the number of days per interval on the smaller timescale
        
        5) num_days_per_larger_interval:

            (int) - the number of days per interval on the larger timescale
        
    Outputs:

        1) patient_count_list:

            (list of 1D numpy arrays) -  a list of all the seizure diaries of random length, each diary can be converted
                                         to a larger time scale

    '''

    # check whether or not the larger intervals can be broken up into discrete smaller intervals
    intervals_are_compatible  = num_days_per_larger_interval % num_days_per_smaller_interval == 0

    # assuming that's true...
    if(intervals_are_compatible == True):

        # calculate the conversion factor between larger intervals and smaller intervals
        interval_conversion_factor = num_days_per_larger_interval/num_days_per_smaller_interval

        # initialize the list of patients with variable seizure diary lengths
        patient_count_list = []

        # for each patient
        for i in range(num_patients):

            # generate the diary length for this specific patient
            num_intervals = np.int_(interval_conversion_factor*np.random.randint(min_num_larger_intervals, max_num_larger_intervals))
            
            # actually generate the patient's seizure diary
            patient_counts = generate_patient_counts(num_intervals, num_days_per_smaller_interval)

            # store the patient's diary
            patient_count_list.append(patient_counts)
    
        return patient_count_list
    
    # otherwise, throw an error
    else:

        raise ValueError('')


def convert_patient_to_larger_interval_timescale(patient_counts, original_num_days_per_interval, new_num_days_per_interval):

    '''

    Sums up the seizure counts from a patient diary such that an array of counts on one timescale are converted to a smaller array
    of counts on a larger timescale (e.g., a diary consisting of 14 daily counts is converted to a diary consisting 2 weekly counts)

    Inputs:

        1) patient_counts:

            (1D numpy array) - the patient seizure diary to be converted

        2) original_num_days_per_interval:

            (int) - the original number of days each interval in the original timescale contained

        3) num_days_per_new_interval:

            (int) - the number of days each interval in the new timescale that the seizure diary will be converted to

    Output:

        1) patient_counts:

            (1D numpy array) - the new seizure diary, which is just the old seizure diary converted to the new timescale

    '''

    # obtain the number of original smaller intervals in the patient diary
    (num_original_intervals,) = patient_counts.shape

    # check if the smaller original intervals are compatible with the new, larger intervals
    are_interval_timescales_compatible = new_num_days_per_interval % original_num_days_per_interval == 0

    # if so, then do the following:
    if( are_interval_timescales_compatible ):

        # find the conversion factor between the number of days 
        interval_conversion_factor = np.int_(new_num_days_per_interval/original_num_days_per_interval)
        are_original_and_new_intervals_compatible = num_original_intervals % interval_conversion_factor == 0

        if( are_original_and_new_intervals_compatible == True ):

            num_new_intervals = np.int_(num_original_intervals/interval_conversion_factor)
            patient_counts = np.sum( np.reshape(patient_counts, (num_new_intervals,  interval_conversion_factor) ), 1 )

            return patient_counts

        else:

            raise ValueError('The patient count data cannot be converted to the specified time interval')
    
    else:

        raise ValueError('The intervals timescales are not compatible')


def generate_fitted_population_parameters(num_patients, min_num_months, max_num_months):
    '''
    
    Generates a population of synthetic patients for the purpose of testing the fitted population parameters

    Inputs:
        1) num_patients:

            (int) - the number of patients needed for the parameter fitting

        2) min_num_months:

            (int) - the minimum number of months that each synthetic patient should have

        3) max_num_months:

            (int) - the maximum number of months that each synthetic patient should have

    Outputs:

        1) median_average_monthly_seizure_frequency: 
        
            (float) - the median average monthly seizure frequency from the generated population of synthetic patients
        
        2) log_log_slope:

            (float) - the slope of the base 10 logarithm of the mean of the bi-weekly seizure count vs the base 10
                      logarithm of the bi-weekly seizure count standard deviations

    '''

    num_days_per_two_weeks = 14
    num_days_per_month = 28

    biweekly_patient_count_list = generate_unrestricted_patient_population(num_patients, min_num_months, max_num_months, num_days_per_two_weeks, num_days_per_month)
    average_monthly_seizure_frequencies = np.zeros(num_patients)
    biweekly_seizure_frequency_log_means = np.zeros(num_patients)
    biweekly_seizure_frequency_log_standard_deviations = np.zeros(num_patients)

    for i in range(num_patients):

        biweekly_patient_counts = biweekly_patient_count_list[i]

        monthly_patient_counts = convert_patient_to_larger_interval_timescale(biweekly_patient_counts, num_days_per_two_weeks, num_days_per_month)
        average_monthly_seizure_frequency = np.mean(monthly_patient_counts)
        average_monthly_seizure_frequencies[i] = average_monthly_seizure_frequency

        biweekly_seizure_frequency_mean = np.mean(biweekly_patient_counts)
        biweekly_seizure_frequency_standard_deviation = np.std(biweekly_patient_counts)
        biweekly_seizure_frequency_log_means[i] = np.log10(biweekly_seizure_frequency_mean)
        biweekly_seizure_frequency_log_standard_deviations[i] = np.log10(biweekly_seizure_frequency_standard_deviation)


    median_average_monthly_seizure_frequency = np.median(average_monthly_seizure_frequencies)
    [log_log_slope, _, _, _, _] = stats.linregress(biweekly_seizure_frequency_log_means, biweekly_seizure_frequency_log_standard_deviations)

    return [median_average_monthly_seizure_frequency, log_log_slope]


def generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect):
    '''

    Generate one patient to fit in an RCT trial

    Inputs:

        1) num_base_weeks:

            (int) - the number of baseline period weeks in the patient's seizure diary

        2) num_test_weeks:

            (int) - the number of testing period weeks in the patient's seizure diary

        3) num_total_base_szs:

            (int) - the minimum number of required seizures throughout baseline 

        4) drug_effect:

            (float) - the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

    Outputs:

        1) weekly_patient_counts

            (1D numpy array) - a seizure diary consisting of baseline and testing counts, with the drug effect applied

    '''

    if(drug_effect < 1):

        num_days_per_week = 7
        num_total_weeks = num_base_weeks + num_test_weeks
        acceptable_baseline = False
        while(not acceptable_baseline):

            weekly_patient_counts = generate_patient_counts(num_total_weeks, num_days_per_week)
            base_weekly_patient_counts = weekly_patient_counts[0:num_base_weeks]
            test_weekly_patient_counts = weekly_patient_counts[num_base_weeks:]

            acceptable_baseline = np.sum(base_weekly_patient_counts) >= num_total_base_szs

        
        for test_week in range(num_test_weeks):

            num_removed = 0
            weekly_count = np.int_(test_weekly_patient_counts[test_week])

            for seizure in range(weekly_count):
                
                if( np.random.uniform() < drug_effect ):
                    
                    num_removed = num_removed + 1

            test_weekly_patient_counts[test_week] = weekly_count  - num_removed

        weekly_patient_counts[num_base_weeks:] = test_weekly_patient_counts

        return weekly_patient_counts
    
    else:

        raise ValueError('The drug effect should be a percentage')


def generate_RCT(num_patients_per_arm, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect, placebo_effect):
    '''

    Generates one RCT with accompanying responses (RR50 and MPC)

    Inputs:
    
        1) num_patients_per_arm:

            (int) - the number of patients per treatment arm

        2) num_base_weeks:

            (int) - the number of baseline period weeks in the patient's seizure diary

        3) num_test_weeks: 

            (int) - the number of testing period weeks in the patient's seizure diary

        4) num_total_base_szs:

            (int) - the minimum number of required seizures throughout baseline 

        5) drug_effect:

            (float) - the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        6)placebo_effect:

            (float) - the percent chance that the ''placebo effect'' will remove one seizure from a patient, applied over all seizures

    Outputs:

        1) placebo_arm_RR50:

            (float) - the 50% responder rate for the placebo arm of the RCT

        2) drug_arm_RR50:

            (float) - the 50% responder rate for the drug arm of the RCT

        3) placebo_arm_MPC:

            (float) - the median percent change for the placebo arm of the RCT

        4) drug_arm_MPC:

            (float) - the median percent change for the drug arm of the RCT

    '''

    placebo_arm_patients = np.zeros([num_patients_per_arm, num_base_weeks + num_test_weeks])
    drug_arm_patients = np.zeros([num_patients_per_arm, num_base_weeks + num_test_weeks])

    for i in range(num_patients_per_arm):

        placebo_arm_patient = generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, placebo_effect)
        drug_arm_patient = generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect + placebo_effect)

        placebo_arm_patients[i, :] = placebo_arm_patient
        drug_arm_patients[i, :] = drug_arm_patient

    baseline_placebo_arm_average_weekly_counts = np.mean(placebo_arm_patients[:, 0:num_base_weeks], 1)
    testing_placebo_arm_average_weekly_counts = np.mean(placebo_arm_patients[:, num_base_weeks:], 1)
    baseline_drug_arm_average_weekly_counts = np.mean(drug_arm_patients[:, :num_base_weeks], 1)
    testing_drug_arm_average_weekly_counts = np.mean(drug_arm_patients[:, num_base_weeks:], 1)

    placebo_arm_percent_change = np.divide(baseline_placebo_arm_average_weekly_counts - testing_placebo_arm_average_weekly_counts, baseline_placebo_arm_average_weekly_counts)
    drug_arm_percent_change = np.divide(baseline_drug_arm_average_weekly_counts - testing_drug_arm_average_weekly_counts, baseline_drug_arm_average_weekly_counts)

    placebo_arm_RR50 = np.sum(placebo_arm_percent_change >= 0.5)/num_patients_per_arm
    drug_arm_RR50 = np.sum(drug_arm_percent_change >= 0.5)/num_patients_per_arm

    placebo_arm_MPC = np.median(placebo_arm_percent_change)
    drug_arm_MPC = np.median(drug_arm_percent_change)

    return [placebo_arm_RR50, drug_arm_RR50, placebo_arm_MPC, drug_arm_MPC]


def simulate_RCT_results(num_patients_per_arm, num_trials, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect, placebo_effect):
    '''

    Calculates the average and standard deviation of the RR50 and MPC responses for all treatment arms over multiple trials

    Inputs:

        1) num_patients_per_arm: 

            (int) - the number of patients per treatment arm

        2) num_trials:

            (int) - the number of trials to calculate the average and standard of the RR50/MPC responses over

        3) num_base_weeks:

            (int) - the number of weeks in the baseline period over all trials

        4) num_test_weeks:

            (int) - the number of weeks in the testing period over all trials

        5) num_total_base_szs: 

            (int) - the minimum required amount of seizures in the baseliner period

        6) drug_effect:

            (float) - the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        7) placebo_effect:

            (float) - the percent chance that the ''placebo effect'' will remove one seizure from a patient, applied over all seizures

    Outputs:

        Technically none, but as the description says, this function prints a string to the console which contains the 
        average and standard deviation of the placebo and drug responses for both RR50 and MPC over multiple trials

    '''

    placebo_arm_RR50_array = np.zeros(num_trials)
    drug_arm_RR50_array = np.zeros(num_trials)
    placebo_arm_MPC_array = np.zeros(num_trials)
    drug_arm_MPC_array = np.zeros(num_trials)

    for i in range(num_trials):

        [placebo_arm_RR50, drug_arm_RR50, placebo_arm_MPC, drug_arm_MPC] = generate_RCT(num_patients_per_arm, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect, placebo_effect)
        placebo_arm_RR50_array[i] = placebo_arm_RR50
        drug_arm_RR50_array[i] = drug_arm_RR50
        placebo_arm_MPC_array[i] = placebo_arm_MPC
        drug_arm_MPC_array[i] = drug_arm_MPC

    placebo_arm_RR50_mean = np.mean(placebo_arm_RR50_array)
    placebo_arm_RR50_std = np.std(placebo_arm_RR50_array)
    drug_arm_RR50_mean = np.mean(drug_arm_RR50_array)
    drug_arm_RR50_std = np.std(drug_arm_RR50_array)
    placebo_arm_MPC_mean = np.mean(placebo_arm_MPC_array)
    placebo_arm_MPC_std = np.std(placebo_arm_MPC_array)
    drug_arm_MPC_mean = np.mean(drug_arm_MPC_array)
    drug_arm_MPC_std = np.std(drug_arm_MPC_array)

    print( '\n\nplacebo RR50: ' + str(np.round(100*placebo_arm_RR50_mean, 2)) + ' +- ' +  str(np.round(100*placebo_arm_RR50_std, 2)) + ' % ' + 
           '\n\ndrug RR50: '    + str(np.round(100*drug_arm_RR50_mean, 2))    + ' +- ' +  str(np.round(100*drug_arm_RR50_std, 2))    + ' % ' +  
           '\n\nplacebo MPC: ' + str(np.round(100*placebo_arm_MPC_mean, 2))  + ' +- ' +  str(np.round(100*placebo_arm_MPC_std, 2))  + ' % ' +  
           '\n\ndrug MPC: '     + str(np.round(100*drug_arm_MPC_mean, 2))     + ' +- ' +  str(np.round(100*drug_arm_MPC_std, 2))     + ' % \n\n' )


num_patients_per_arm = 153
num_trials = 100
num_base_weeks = 8
num_test_weeks = 12
num_total_base_szs = 4
drug_effect = 0.2
placebo_effect = 0

simulate_RCT_results(num_patients_per_arm, num_trials, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect, placebo_effect)

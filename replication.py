import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


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
    '''
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

        # check if the time length of the sieuzre diary allows for conversion to the larger timescale
        are_original_and_new_intervals_compatible = num_original_intervals % interval_conversion_factor == 0

        # if the seizure diary length does allow for sucha thing. then...
        if( are_original_and_new_intervals_compatible == True ):

            # calcualte the number of larget intervals in the new timescale for an equivalent seziure diary
            num_new_intervals = np.int_(num_original_intervals/interval_conversion_factor)

            # actually create the summarized seizure diary on the larger time-scale
            patient_counts = np.sum( np.reshape(patient_counts, (num_new_intervals,  interval_conversion_factor) ), 1 )

            return patient_counts

        # ifthe seizure diary length is, in fact, not compatible, then...
        else:

            # send an error message to the user saying that the patient count data isn't compatible
            raise ValueError('The patient count data cannot be converted to the specified time interval')
    
    # if the timescale are not compatible, then...
    else:

        # send a message to the user saying that the time-scales aren't compatible
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

    # initiailize the timescales needed for calculating the bi-weekly log-log slope and monthly median freq.
    num_days_per_two_weeks = 14
    num_days_per_month = 28

    # Generate all the patients needed for one estimate of the log-log slope and one estimate of the median frequency.
    # These patients should have seizure diaries with counts on a bi-weekly time interval and a diary length that allows
    # for conversion to a monthly time scale.
    biweekly_patient_count_list = generate_unrestricted_patient_population(num_patients, min_num_months, max_num_months, num_days_per_two_weeks, num_days_per_month)
    
    # initialize the arrays that will store the information needed to calculate the slope and median
    average_monthly_seizure_frequencies = np.zeros(num_patients)
    biweekly_seizure_frequency_log_means = np.zeros(num_patients)
    biweekly_seizure_frequency_log_standard_deviations = np.zeros(num_patients)

    # for every patient generated...
    for i in range(num_patients):

        # take their seziure diary
        biweekly_patient_counts = biweekly_patient_count_list[i]

        # convert their bi-weekly seizure diary to a monthly time scale
        monthly_patient_counts = convert_patient_to_larger_interval_timescale(biweekly_patient_counts, num_days_per_two_weeks, num_days_per_month)

        # calculate and store the monthly seizure frequency
        average_monthly_seizure_frequency = np.mean(monthly_patient_counts)
        average_monthly_seizure_frequencies[i] = average_monthly_seizure_frequency

        # calculate the base 10 logarithms of the bi-weekly mean and bi-weekly standard deviations
        biweekly_seizure_frequency_mean = np.mean(biweekly_patient_counts)
        biweekly_seizure_frequency_standard_deviation = np.std(biweekly_patient_counts)
        biweekly_seizure_frequency_log_means[i] = np.log10(biweekly_seizure_frequency_mean)
        biweekly_seizure_frequency_log_standard_deviations[i] = np.log10(biweekly_seizure_frequency_standard_deviation)

    # find the median monhtly seizure frequency amongst all the generated patients
    median_average_monthly_seizure_frequency = np.median(average_monthly_seizure_frequencies)

    # use linear regression to find the slope of the line between the base-10 logarithm of the bi-weekly mean
    # and the base-10 logarithm of the bi-weekly standard deviation
    [log_log_slope, _, _, _, _] = stats.linregress(biweekly_seizure_frequency_log_means, biweekly_seizure_frequency_log_standard_deviations)

    return [median_average_monthly_seizure_frequency, log_log_slope]


def generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu, drug_effect_sigma):
    '''

    Generate one patient to fit in an RCT trial

    Inputs:

        1) num_base_weeks:

            (int) - the number of baseline period weeks in the patient's seizure diary

        2) num_test_weeks:

            (int) - the number of testing period weeks in the patient's seizure diary

        3) num_total_base_szs:

            (int) - the minimum number of required seizures throughout baseline 

        4) drug_effect_mu:

            (float) - the mean of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        5) drug_effect_sigma:

            (float) - the standard deviation of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

    Outputs:

        1) weekly_patient_counts

            (1D numpy array) - a seizure diary consisting of baseline and testing counts, with the drug effect applied

    '''

    # check if the 95 % confidence interval of the drug effect parameters are outside of the interval [-100 % , 100 %]
    if( ( ( drug_effect_mu + 1.96*drug_effect_sigma ) < 1 ) and ( ( drug_effect_mu - 1.96*drug_effect_sigma ) > -1 ) ):

        # if not, then randomly generate the drug effect according to a gaussian distribution
        '''

        There's some kind of problem going with the drug effect: using the mean gives perfectly fine results, but using the gausssian generated drug effect gives strange results

        '''
        drug_effect = np.random.normal(drug_effect_mu, drug_effect_sigma)

        # if the drug effect is a percentage over 100%, then...
        if( drug_effect <= 1 ):

            # change the drug_effect to be below 100 %
            drug_effect = 0.999999999

        # if the drug effect is a percentage below -100%, then...
        if( drug_effect >= -1 ):

            # change the drug_effect to be above -100 %
            drug_effect = -0.999999999

        # initialize the information needed to generate a weekly seizure diary 
        # with a specified baseline period and testing period
        num_days_per_week = 7
        num_total_weeks = num_base_weeks + num_test_weeks
        
        # use a while loop to make sure that any generated patient has the minimum required baseline seizures
        # if they don't, discard their data and make a new diary
        # if they do, store their data and move on
        acceptable_baseline = False
        while(not acceptable_baseline):

            weekly_patient_counts = generate_patient_counts(num_total_weeks, num_days_per_week)
            base_weekly_patient_counts = weekly_patient_counts[0:num_base_weeks]
            test_weekly_patient_counts = weekly_patient_counts[num_base_weeks:]

            acceptable_baseline = np.sum(base_weekly_patient_counts) >= num_total_base_szs

        # for each week in the testing period:
        for test_week in range(num_test_weeks):

            # initialize the number of seizures removed for this week in the testing period
            num_removed = 0

            # get the relevant weekly count from the testing period
            weekly_count = np.int_(test_weekly_patient_counts[test_week])

            # for each seizure in the relevant weekly seizure counts...
            for seizure in range(weekly_count):
                
                # if the drug effect is positive:
                if(drug_effect >= 0):

                    # generate a percentage between 0 % and 100 %
                    # if that percentage is less than the drug effect size:
                    if( np.random.uniform(0, 1) <= drug_effect ):
                    
                        # say that a seizure has been removed
                        num_removed = num_removed + 1

                # if the drug effect is negative...
                if(drug_effect < 0):

                    # generate a percentage between 0 % and -100 %
                    # if that percentage is greater than the drug effect size...
                    if( np.random.uniform(-1, 0) >= drug_effect ):

                        # say that a seizure has been added
                        num_removed = num_removed - 1

            # actually remove/add the number of seizures as determined by the random percentages
            test_weekly_patient_counts[test_week] = weekly_count  - num_removed

        # replace the old weekly testing period count with the new drug-modified weekly seizure count
        weekly_patient_counts[num_base_weeks:] = test_weekly_patient_counts

        return weekly_patient_counts
    
    # else if the drug effect is not a percentage:
    else:

        # let the user know that their drug effect value doesn't make sense
        raise ValueError('The 95% confidence interval of the drug effect parameters should not be outside the interval [-100 % , 100 %]')


def generate_RCT(num_patients_per_arm, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu, drug_effect_sigma, placebo_effect_mu):
    '''

    Generates one RCT (1 drug arm and 1 placebo arm) with accompanying responses (RR50 and MPC)

    Inputs:
    
        1) num_patients_per_arm:

            (int) - the number of patients per treatment arm

        2) num_base_weeks:

            (int) - the number of baseline period weeks in the patient's seizure diary

        3) num_test_weeks: 

            (int) - the number of testing period weeks in the patient's seizure diary

        4) num_total_base_szs:

            (int) - the minimum number of required seizures throughout baseline 

        5) drug_effect_mu:

            (float) - the mean of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        6) drug_effect_sigma:

            (float) - the standard deviation of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        7) placebo_effect_mu:

            (float) - the mean of the percent chance that the 'placebo effect will remove one seizure from a patient, applied over all seizures

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

    # initialize the 2D numpy arrays that will store the weekly seizure diaries for all the 
    # patients in the drug arm and placebo arm, with 1 2D numpy array per arm
    placebo_arm_patients = np.zeros([num_patients_per_arm, num_base_weeks + num_test_weeks])
    drug_arm_patients = np.zeros([num_patients_per_arm, num_base_weeks + num_test_weeks])

    # for each specified patient:
    for i in range(num_patients_per_arm):

        # generate one patient for both arms
        placebo_arm_patient = generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, placebo_effect_mu, drug_effect_sigma)
        drug_arm_patient = generate_RCT_patient(num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu + placebo_effect_mu, drug_effect_sigma)

        # store the generated synthetic patients into their respective arms
        placebo_arm_patients[i, :] = placebo_arm_patient
        drug_arm_patients[i, :] = drug_arm_patient

    # for each arm, separate each patients weekly counts into a baseline and testing period
    baseline_placebo_arm_average_weekly_counts = np.mean(placebo_arm_patients[:, 0:num_base_weeks], 1)
    testing_placebo_arm_average_weekly_counts = np.mean(placebo_arm_patients[:, num_base_weeks:], 1)
    baseline_drug_arm_average_weekly_counts = np.mean(drug_arm_patients[:, 0:num_base_weeks], 1)
    testing_drug_arm_average_weekly_counts = np.mean(drug_arm_patients[:, num_base_weeks:], 1)

    # avoid divide-by-zero error by making zero-valued baseline frequencies very small non-zero values
    for i in range(num_patients_per_arm):
        if( baseline_placebo_arm_average_weekly_counts[i] == 0):
            baseline_placebo_arm_average_weekly_counts[i] = 0.00000001
        if( baseline_drug_arm_average_weekly_counts[i] == 0):
            baseline_drug_arm_average_weekly_counts[i] = 0.00000001

    # calculate the percent change for each patient in both arms
    placebo_arm_percent_change = np.divide(baseline_placebo_arm_average_weekly_counts - testing_placebo_arm_average_weekly_counts, baseline_placebo_arm_average_weekly_counts)
    drug_arm_percent_change = np.divide(baseline_drug_arm_average_weekly_counts - testing_drug_arm_average_weekly_counts, baseline_drug_arm_average_weekly_counts)

    # calculate the 50% responder rate over all the percent changes in the placebo arm as well as the drug arm
    placebo_arm_RR50 = np.sum(placebo_arm_percent_change >= 0.5)/num_patients_per_arm
    drug_arm_RR50 = np.sum(drug_arm_percent_change >= 0.5)/num_patients_per_arm

    # calculate the median percent change over all the percent changes in the placebo arm as well as the drug arm
    placebo_arm_MPC = np.median(placebo_arm_percent_change)
    drug_arm_MPC = np.median(drug_arm_percent_change)

    return [placebo_arm_RR50, drug_arm_RR50, placebo_arm_MPC, drug_arm_MPC]


def simulate_RCT_results(num_patients_per_arm, num_trials, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu, drug_effect_sigma, placebo_effect_mu):
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

        5) drug_effect_mu:

            (float) - the mean of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        7) drug_effect_sigma:

            (float) - the standard deviation of the percent chance that the simulated drug will remove one seizure from a patient, applied over all seizures

        8) placebo_effect_mu:

            (float) - the mean of the percent chance that the 'placebo effect will remove one seizure from a patient, applied over all seizures

    Outputs:

        Technically none, but as the description says, this function prints a string to the console which contains the 
        average and standard deviation of the placebo and drug responses for both RR50 and MPC over multiple trials

    '''

    # initialize the arrays that will be used to store the 50% responder rate and median percent changes from all the specified trials
    placebo_arm_RR50_array = np.zeros(num_trials)
    drug_arm_RR50_array = np.zeros(num_trials)
    placebo_arm_MPC_array = np.zeros(num_trials)
    drug_arm_MPC_array = np.zeros(num_trials)

    # over all the specified trials:
    for i in range(num_trials):

        # generate the 50% responder rate and median percent change for both the placebo arm and drug arm of one trial
        [placebo_arm_RR50, drug_arm_RR50, placebo_arm_MPC, drug_arm_MPC] = \
                generate_RCT(num_patients_per_arm, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu, drug_effect_sigma, placebo_effect_mu)
        
        # store the aformentioned values
        placebo_arm_RR50_array[i] = placebo_arm_RR50
        drug_arm_RR50_array[i] = drug_arm_RR50
        placebo_arm_MPC_array[i] = placebo_arm_MPC
        drug_arm_MPC_array[i] = drug_arm_MPC

    # calculate the mean and standard deviation of the 50% responder rate and median percent change for the drug and placebo arm over all the trials
    placebo_arm_RR50_mean = np.mean(placebo_arm_RR50_array)
    placebo_arm_RR50_std = np.std(placebo_arm_RR50_array)
    drug_arm_RR50_mean = np.mean(drug_arm_RR50_array)
    drug_arm_RR50_std = np.std(drug_arm_RR50_array)
    placebo_arm_MPC_mean = np.mean(placebo_arm_MPC_array)
    placebo_arm_MPC_std = np.std(placebo_arm_MPC_array)
    drug_arm_MPC_mean = np.mean(drug_arm_MPC_array)
    drug_arm_MPC_std = np.std(drug_arm_MPC_array)

    return [placebo_arm_RR50_mean, placebo_arm_RR50_std, drug_arm_RR50_mean, drug_arm_RR50_std, 
            placebo_arm_MPC_mean,  placebo_arm_MPC_std,  drug_arm_MPC_mean,  drug_arm_MPC_std]

# log-log plot and histogram parameters
num_patients = 10000
min_num_months = 30
max_num_months = 128
decimal_round = 2
[median_average_monthly_seizure_frequency, log_log_slope] = generate_fitted_population_parameters(num_patients, min_num_months, max_num_months)

print(          '\n\nlog-log slope: '          +                str(np.round(log_log_slope, decimal_round))             + 
      '\n\nmedian monthly seizure frequency: ' + str(np.round(median_average_monthly_seizure_frequency, decimal_round)) + '\n\n')

# parameters for calculating the mean and standard deviation of the 
# 50% responder rate and median percent change for both the drug and
# placebo arms over many RCTs
num_patients_per_arm = 153
num_trials = 5000
num_base_weeks = 8
num_test_weeks = 12
num_total_base_szs = 4
drug_effect_mu = 0.2
drug_effect_sigma = 0.05
placebo_effect_mu = 0

[placebo_arm_RR50_mean, placebo_arm_RR50_std, drug_arm_RR50_mean, drug_arm_RR50_std, 
 placebo_arm_MPC_mean,  placebo_arm_MPC_std,  drug_arm_MPC_mean,  drug_arm_MPC_std] = \
            simulate_RCT_results(num_patients_per_arm, num_trials, num_base_weeks, num_test_weeks, num_total_base_szs, drug_effect_mu, drug_effect_sigma, placebo_effect_mu)

print( '\n\nplacebo RR50: ' + str(np.round(100*placebo_arm_RR50_mean, decimal_round)) + ' +- ' +  str(np.round(100*placebo_arm_RR50_std, decimal_round)) + ' % ' + 
       '\n\ndrug RR50: '    + str(np.round(100*drug_arm_RR50_mean, decimal_round))    + ' +- ' +  str(np.round(100*drug_arm_RR50_std, decimal_round))    + ' % ' +  
       '\n\nplacebo MPC: '  + str(np.round(100*placebo_arm_MPC_mean, decimal_round))  + ' +- ' +  str(np.round(100*placebo_arm_MPC_std, decimal_round))  + ' % ' +  
       '\n\ndrug MPC: '     + str(np.round(100*drug_arm_MPC_mean, decimal_round))     + ' +- ' +  str(np.round(100*drug_arm_MPC_std, decimal_round))     + ' % \n\n' )

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import textwrap
import os
import json
import time


def generate_patient(shape, scale, alpha, beta,
                     time_scale_conversion, 
                     num_baseline_intervals, num_total_intervals, 
                     min_required_baseline_seizure_count):
    '''

    This function randomly generate one sythetic patient's seizure diary as according to the NV model.

    Inputs:

        1) shape:

            (float) - first group level parameter for the NV model
        
        2) scale:

            (float) - second group level parameter for the NV model
        
        3) alpha:

            (float) - third group level parameter for the NV model
        
        4) beta:

            (float) - fourth group level parameter for the NV model
        
        5) time_scale_conversion:
        
            (float) -  A time-scale conversion factor which determines what period of time the 
                       
                       seizure counts will be generated in (i.e., hourly, daily, weekly, biweeekly, monthly, etc.).

                       The group level parameters are optimized for daily seizure counts, so in order to generate daily

                       seizure counts, just set this parameter equal. To generate weekly counts, set this parameter equal to

                       7, and to generate hourly counts, set this parameter equal to 1/24. In general, to convert from a daily 
                       
                       timescale to a larger time scale, set this parameter equal to the number of days in each interval of that 

                       time scale, and to convert to a smaller time scale, set this parameter equal to 1 divided by the number of

                       intervals in 1 day.

        6) num_baseline_intervals:

            (int) - the number of intervals in the baseline period of each patient
        
        7)  num_total_intervals:

            (int) - the number of intervals in each patient's seizure diary
        
        8) min_required_baseline_seizure_count:

            (int) - eligibility criteria which is the minimum amount of seizures that have to be in the baseline of each patient
    
    Outputs:

        1) patient_seizure_counts:

            (1D NUmpy array) - One patient's seizure diary, generated as according to the NV model

    '''
    
    # initialize the boolean status of the baseline rate as not acceptable according to eligibility criteria
    acceptable_baseline_rate = False
    
    # initialize the array which will contain the seizure diary for one patient
    patient_seizure_counts = np.zeros(num_total_intervals)
    
    # while the baseline rate is not not acceptable according to eligibility criteria:
    while(not acceptable_baseline_rate):
    
        # randomly generate n and p for each different patient
        n = np.random.gamma(shape, 1/scale)
        p = np.random.beta(alpha, beta)
    
        # convert the n and p parameters to mean and overdispersion parameters 
        mean = n*( (1 - p)/p )
        overdispersion = 1/n
    
        # for each interval in the seizure diary:
        for interval_index in range(num_total_intervals):
        
            # randomly generate seizure counts as according to gamma-poisson mixture (continuous equivalent of negative binomial)
            rate = np.random.gamma(time_scale_conversion/overdispersion, mean*overdispersion)
            count = np.random.poisson(rate)
        
            # store the seizure count in the patient's seizure diary
            patient_seizure_counts[interval_index] = count
        
        # if the number of intervals in the baseline period is not zero, and the baseline rate satisfies the eligibility criteria, then:
        if( (num_baseline_intervals != 0) and ( np.sum(patient_seizure_counts[0:num_baseline_intervals]) >= min_required_baseline_seizure_count ) ):
            
            # say that the baseline rate is acceptable
            acceptable_baseline_rate = True
        
        # else if the baseline rate is zero:
        elif( num_baseline_intervals == 0 ):
        
            # say that the baseline rate is acceptable, since it doesn't really matter
            acceptable_baseline_rate = True
    
    return patient_seizure_counts


def get_patient_statistical_features(shape, scale, alpha, beta,
                                     min_num_weeks, max_num_weeks):
    '''

    This function generates one synthetic patient's weekly seizure diary as according to the NV model, and calculates 
    
    their statistical features, as described in the outputs. This patient's seizure diary length is randomly determined

    by a discrete uniform distribution whose minimum and maximum number of weeks is determined by the user.

    Inputs:

        1) shape:
        
            (float) - first group level parameter for the NV model

        2) scale:

            (float) - second group level parameter for the NV model

        3) alpha:

            (float) - third group level parameter for the NV model

        4) beta:

            (float) - fourth group level parameter for the NV model

        5) min_num_weeks:

            (int) - the minimum number of weeks that any generated synthetic patient can have

        6) max_num_weeks:

            (int) - the maximum number of weeks that any generated synthetic patient can have

    Outputs:

        1) biweekly_log10mean:

            (float) - the base 10 logarithm of the biweekly seizure count mean of one synthetic patient
        
        2) biweekly_log10std:

            (float) - the base 10 logarithm of the biweekly seizure count standard deviation of one synthetic patient
        
        3) monthly_seizure_frequency:

            (float) - the monthly seizure count mean of one syntheic patient

    '''

    # calculate the length of one synthetic patient's weekly seizure diary
    num_testing_intervals = np.random.randint(min_num_weeks, max_num_weeks)

    # generate one synthetic patient's weekly seizure diary according to NV model and the randomly generated diary length
    patient_weekly_seizure_counts = generate_patient(shape, scale, alpha, beta,
                                                      7, 0, num_testing_intervals, 0)
    
    # convert the weekly seizure diary into a biweekly seizure diary
    twoWeeks = []
    for y in range(0, len(patient_weekly_seizure_counts)-1, 2):
        twoWeeks.append(patient_weekly_seizure_counts[y]+patient_weekly_seizure_counts[y+1])
    patient_biweekly_seizure_counts = np.array(twoWeeks)

    # convert the weekly seizure diary into a monthly seizure diary
    months = []
    for y in range(0, len(patient_weekly_seizure_counts)-3, 4):
        months.append(patient_weekly_seizure_counts[y]+patient_weekly_seizure_counts[y+1]+patient_weekly_seizure_counts[y+2]+patient_weekly_seizure_counts[y+3])
    patient_monthly_seizure_counts = np.array(months)

    # calculate the base 10 logarithm of the biweekly seizure count mean
    biweekly_log10mean = np.log10(np.mean(patient_biweekly_seizure_counts))

    # calculate the base 10 logarithm of the biweekly seizure count mean
    biweekly_log10std = np.log10(np.std(patient_biweekly_seizure_counts))

    # calculate the omnthly seizure count mean
    monthly_seizure_frequency = np.mean(patient_monthly_seizure_counts)
    
    return [biweekly_log10mean, biweekly_log10std, monthly_seizure_frequency]


def get_patient_population_statistical_features(shape, scale, alpha, beta, 
                                                min_num_weeks, max_num_weeks, 
                                                num_patients):
    '''

    This function generates one synthetic patient population. The patients have weekly seizure diaries which

    are not all of equal length: rather, each patient's weekly seizure diary length is randomly generated 

    according to a uniform distribution, and the minimum as well as the maximum number of weeks in this

    uniform distribution is specified by the user.

    Inputs:

        1) shape:
        
            (float) - first group level parameter for the NV model
        
        2) scale:

            (float) - second group level parameter for the NV model

        3) alpha:

            (float) - third group level parameter for the NV model
        
        4) beta:

            (float) - fourth group level parameter for the NV model

        5) min_num_weeks:

            (int) - the minimum number of weeks that any patient in the synthetic patient population can have
        
        6) max_num_weeks:

            (int) -  the maximum number of weeks that any patient in the synthetic patient population can have
        
        7) num_patients:

            (int) - the total number of patients in the synthetic patient population
    
    Outputs:

        1) median_monthly_seizure_frequency:
        
            (float) - the median monthly seizure frequency as estimated from one generated synthetic patient population

        2) log_log_slope:

            (float) - the slope of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
            
                      seizure count standard deviation

        3) log_log_intercept:

            (float) - the intercept of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation

        4) r_value:

            (float) - the correlation of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation

        5) monthly_seizure_frequencies:

            (1D Numpy array) - the array of mnthly seizure frequencies from which the median monthly seizure frequency 
                               
                               was estimated
        
        6) biweekly_log10means:
        
            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count means used to construct 
                               
                               the log-log plot

        biweekly_log10stds

            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count standard deviations used to 
            
                               construct the log-log plot

    '''

    # initialize the array of base 10 logarithms of biweekly seizure count means
    biweekly_log10means = np.zeros(num_patients)

    # initialize the array of base 10 logarithms of biweekly seizure count standard deviations
    biweekly_log10stds = np.zeros(num_patients)

    # initialize the array of monthly seizure frequencies
    monthly_seizure_frequencies = np.zeros(num_patients)

    # for each patient in the synthetic patient population:
    for patient_index in range(num_patients):
    
        # calculate each individual patient's statistical features:
        #   1) base 10 logarithm of biweekly mean
        #   2) base 10 logarithm of biweekly standard deviation
        #   3) monthly seizure count mean
        [biweekly_log10mean, biweekly_log10std, monthly_seizure_frequency] = \
                        get_patient_statistical_features(shape, scale, alpha, beta,
                                                         min_num_weeks, max_num_weeks)
    
        # store the patient's statistical features
        biweekly_log10means[patient_index] = biweekly_log10mean
        biweekly_log10stds[patient_index] = biweekly_log10std
        monthly_seizure_frequencies[patient_index] = monthly_seizure_frequency

    # calculate the median monthly seizure frequency
    median_monthly_seizure_frequency = np.median(monthly_seizure_frequencies)

    # use linear regression to estimate the slope, intercept and correlation coefficient of log-log plot of biweekly mean vs biweekly standard deviation
    [log_log_slope, log_log_intercept, r_value, _, _] = stats.linregress(biweekly_log10means, biweekly_log10stds)
    
    return [median_monthly_seizure_frequency, log_log_slope, log_log_intercept, r_value, monthly_seizure_frequencies, biweekly_log10means, biweekly_log10stds]


def apply_effect(effect_mu, effect_sigma,
                 patient_seizure_counts,
                 num_baseline_intervals, num_testing_intervals):
    '''

    This function modifies the seziure counts in the testing period of one patient's seizure diary 

    according to a randomly generated effect size (generated via Normal distribution). If the effect is

    postive, it removes seizures. If it is negative, it adds them.

    Inputs:

        1) effect_mu:

            (float) - the mean of the desired effect's distribution
        
        2) effect_sigma:

            (float) - the standard deviation of the desired effect's distribution

        3) patient_seizure_counts:
        
            (1D Numpy array) - the seizure diary of one patient

        4) num_baseline_intervals:
        
            (int) - the number of intervals in the baseline period of each patients

        5) num_testing_intervals:

            (int) - the number of intervals in the testing period of each patient
        
    Outputs:
    
        1) patient_seizure_counts:

            (1D Numpy array) - the seizure diary of one patient, with seizure counts of the testing period modified by the effect

    '''

    # randomly generate an effect from the given distribution
    effect = np.random.normal(effect_mu, effect_sigma)
    
    # if the effect is greater than 100%:
    if(effect > 1):
        
        # threshold it so that it cannot be greater than 100%
        effect = 1
    
    # extract only the seizure counts from the testing period
    testing_counts = patient_seizure_counts[num_baseline_intervals:]
    
    # for each seizure count in the testing period:
    for testing_count_index in range(num_testing_intervals):
        
        # initialize the number of seizures removed per each seizure count
        num_removed = 0

        # extract only the relevant seizure count from the testing period
        testing_count = testing_counts[testing_count_index]
            
        # for each individual seizure in the relevant seizure count
        for count_index in range(np.int_(testing_count)):
                
            # if a randomly generated number between 0 and 1 is less than or equal to the absolute value of the effect
            if(np.random.random() <= np.abs(effect)):
                    
                # iterate the number of seizures removed by the sign of the effect
                num_removed = num_removed + np.sign(effect)
        
        # remove (or add) the number of seizures from the relevant seizure count as determined by the probabilistic algorithm above
        testing_counts[testing_count_index] = testing_count  - num_removed
    
    # set the patient's seizure counts from the testing period equal to the newly modified seizure counts
    patient_seizure_counts[num_baseline_intervals:] = testing_counts
    
    return patient_seizure_counts


def calculate_endpoints(seizure_counts, num_baseline_intervals, num_patients_per_arm):
    
    '''

    This function calculates the 50% responder rate and median percent change for all the 

    patients in the arm of one trial

    Inputs:

        1) seizure_counts:

            (1D Numpy array) - an array of seizure diaries from all patients of one arm of one trial
        
        2) num_baseline_intervals:

            (int) - the number of intervals in the baseline period of each patient
        
        3) num_patients_per_arm:

            (int) - the number of patients for both the placebo arm and the drug arm of a trial

    Outputs:

        1) RR50:
        
            (float) - the 50% responder rate of all the patients in the arm of one trial
        
        2) MPC:

            (float) - the median percent change of all the patients in the arm of one trial

    '''

    # separate the seizure counts into baseline and testing periods
    baseline_seizure_counts = seizure_counts[:, 0:num_baseline_intervals]
    testing_seizure_counts = seizure_counts[:, num_baseline_intervals:]

    # calculate the seizure frequencies of both periods for all patients
    baseline_seizure_frequencies = np.mean(baseline_seizure_counts, 1)
    testing_seizure_frequencies = np.mean(testing_seizure_counts, 1)

    # for each patient in one arm of one trial 
    for patient_index in range(num_patients_per_arm):

        # if that patient's baseline seizure frequency is zero:
        if( baseline_seizure_frequencies[patient_index] == 0 ):

            # set the baselinse seizure frequency for that patient equal to a very small positive number
            # this is done for the purpose of avoiding divide-by-zero errors when calculating percent change
            baseline_seizure_frequencies[patient_index] = 0.00000001

    # calculate the percent change between the baseline and testing periods for all patients
    percent_changes = np.divide(baseline_seizure_frequencies - testing_seizure_frequencies, baseline_seizure_frequencies)
    
    # calculate the 50% responder rate and median percent change across all the percent changes
    RR50 = 100*np.sum(percent_changes >= 0.5)/num_patients_per_arm
    MPC = 100*np.median(percent_changes)
    
    return [RR50, MPC]


def calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                  placebo_effect_mu, placebo_effect_sigma,
                                  drug_effect_mu, drug_effect_sigma, drug_arm_flag,
                                  num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                  time_scale_conversion, min_required_baseline_seizure_count, num_trials):

    '''

    This function calculates the mean and standard deviation of: 

        1) the 50% responder rate of one arm of a trial over N (N = num_trials) trials

        2) the median percent change of one arm of a trial over N (N = num_trials) trials
    
    The trial arm in question can either be the placebo arm or the drug arm, depending on the drug_arm_flag parameter.

    The mean and standard deviation of each endpoint (RR50, MPC) over multiple trial are referred to within this code 
    
    documentation as the statistics of that endpoint.

    Inputs:

        1) shape:

            (float) - first group level parameter for the NV model
        
        2) scale: 
        
            (float) - second group level parameter for the NV model

        3) alpha:

            (float) - third group level parameter for the NV model
        
        4) beta:

            (float) - fourth group level parameter for the NV model
 
        5) placebo_effect_mu:

            (float) - the mean of the normally distributed placebo effect distribution
        
        6) placebo_effect_sigma:

            (float) - the standard deviation of the normally distributed placebo effect distribution

        7) drug_effect_mu:

            (float) - the mean of the normally distributed drug efficacy distribution
        
        8) drug_effect_sigma:

            (float) - the standard deviation of the normally distributed drug efficacy distribution
        
        9) drug_arm_flag:

            (boolean) - a boolean that determeins whether or not the endpoints are being calculated 

                        over the drug arm or the placebo arm

        10) num_patients_per_arm:

            (int) - the number of patients for both the placebo arm and the drug arm of one trial
        
        11) num_baseline_intervals:

            (int) - the number of intervals in the baseline period of each patient
        
        12) num_testing_intervals:

            (int) - the number of intervals in the testing period of each patient
 
        13) time_scale_conversion:
        
            (float) -  A time-scale conversion factor which determines what period of time the 
                       
                       seizure counts will be generated in (i.e., hourly, daily, weekly, biweeekly, monthly, etc.).

                       The group level parameters are optimized for daily seizure counts, so in order to generate daily

                       seizure counts, just set this parameter equal. To generate weekly counts, set this parameter equal to

                       7, and to generate hourly counts, set this parameter equal to 1/24. In general, to convert from a daily 
                       
                       timescale to a larger time scale, set this parameter equal to the number of days in each interval of that 

                       time scale, and to convert to a smaller time scale, set this parameter equal to 1 divided by the number of

                       intervals in 1 day.
        
        14) min_required_baseline_seizure_count:

            (int) - eligibility criteria which is the minimum amount of seizures that have to be in the baseline of each patient
        
        15) num_trials:

            (int) -  the number of trials that the endpoints are going to be averaged over
    
    Outputs:

        1) RR50_mean:

            (float) - the mean of the 50% responder rate of one arm of a trial over N (N = num_trials) trials
        
        2) RR50_std:

            (float) - the standard deviation of the 50% responder rate of one arm of a trial over N (N = num_trials) trials
        
        3) MPC_mean:

            (float) - the mean of the median percent change of one arm of a trial over N (N = num_trials) trials
        
        4) MPC_std:

            (float) - the standard deviation of the median percent change of one arm of a trial over N (N = num_trials) trials

    '''

    # calculate the total number of intervals in each patient's seizure diary
    num_total_intervals = num_baseline_intervals + num_testing_intervals

    # initialize the array of 50% responder rates for each trial
    RR50_array = np.zeros(num_trials)

    # initialize the array of median percent changes for each trial
    MPC_array = np.zeros(num_trials)

    # for each simulated trial:
    for trial_index in range(num_trials):
        
        # initialize the 2D numpy array which contains the seizure diaries of all patients within a trial arm
        seizure_counts = np.zeros((num_patients_per_arm, num_total_intervals))

        # for each patient within the arm of one trial:
        for patient_index in range(num_patients_per_arm):
        
            # generate the seizure diary of one patient
            seizure_counts[patient_index, :] = generate_patient(shape, scale, alpha, beta,
                                                                time_scale_conversion, num_baseline_intervals, 
                                                                num_total_intervals, min_required_baseline_seizure_count)
        
            # for that patient, apply a placebo effect which modifies their seizure counts
            seizure_counts[patient_index, :] = apply_effect(placebo_effect_mu, placebo_effect_sigma,
                                                            seizure_counts[patient_index, :],
                                                            num_baseline_intervals, num_testing_intervals)
            
            # if the trial arm is specified as the drug arm of a trial:
            if(drug_arm_flag):

                # apply a drug efficacy which modifies the patient's seizure counts
                seizure_counts[patient_index, :] = apply_effect(drug_effect_mu, drug_effect_sigma,
                                                                seizure_counts[patient_index, :],
                                                                num_baseline_intervals, num_testing_intervals)

        # calculate the 50% responder rate and median percent change of all the patients in one trial arm
        [RR50, MPC] = calculate_endpoints(seizure_counts, num_baseline_intervals, num_patients_per_arm)
    
        # store the 50% responder rate and median percent change within their respective arrays
        RR50_array[trial_index] = RR50
        MPC_array[trial_index] = MPC

    # calculate the mean and standard deviation of the 50% responder rates over all the trials
    RR50_mean = np.mean(RR50_array)
    RR50_std = np.std(RR50_array)

    # calculate the mean and standard deviation of the median percent changes over all the trials
    MPC_mean = np.mean(MPC_array)
    MPC_std = np.std(MPC_array)

    return [RR50_mean, RR50_std, MPC_mean, MPC_std]


def calculate_placebo_and_drug_arm_endpoint_statistics(shape, scale, alpha, beta, 
                                                       placebo_effect_mu, placebo_effect_sigma,
                                                       drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       time_scale_conversion, min_required_baseline_seizure_count, num_trials):
    '''

    This function calculates the mean and standard deviation of: 

        1) the 50% responder rate of the placebo arm over N (N = num_trials) trials

        2) the median percent change of the placebo arm over N (N = num_trials) trials

        3) the 50% responder rate of the drug arm N (N = num_trials) trials

        4) the median percent change of the drug arm N (N = num_trials) trials
    
    The mean and standard deviation of each endpoint (RR50, MPC) over multiple trial are referred to within this code 
    
    documentation as the statistics of that endpoint.

    Inputs:

        1) shape:

            (float) - first group level parameter for the NV model

        2) scale:
        
            (float) - second group level parameter for the NV model
        
        3) alpha:
            
            (float) - third group level parameter for the NV model
        
        4) beta:

            (float) - fourth group level parameter for the NV model
 
        5) placebo_effect_mu:
        
            (float) - the mean of the normally distributed placebo effect distribution
        
        6) placebo_effect_sigma:

            (float) - the standard deviation of the normally distributed placebo effect distribution
        
        7) drug_effect_mu:
        
            (float) - the mean of the normally distributed drug efficacy distribution
        
        8)  drug_effect_sigma:

            (float) - the standard deviation of the normally distributed drug efficacy distribution

        9) num_patients_per_arm:
        
            (int) - the number of patients for both the placebo arm and the drug arm of one trial
        
        10) num_baseline_intervals:
        
            (int) - the number of intervals in the baseline period of each patient
        
        11) num_testing_intervals:

            (int) - the number of intervals in the testing period of each patient

        12) time_scale_conversion:
        
            (float) -  A time-scale conversion factor which determines what period of time the 
                       
                       seizure counts will be generated in (i.e., hourly, daily, weekly, biweeekly, monthly, etc.).

                       The group level parameters are optimized for daily seizure counts, so in order to generate daily

                       seizure counts, just set this parameter equal. To generate weekly counts, set this parameter equal to

                       7, and to generate hourly counts, set this parameter equal to 1/24. In general, to convert from a daily 
                       
                       timescale to a larger time scale, set this parameter equal to the number of days in each interval of that 

                       time scale, and to convert to a smaller time scale, set this parameter equal to 1 divided by the number of

                       intervals in 1 day.

        13) min_required_baseline_seizure_count:
        
            (int) - eligibility criteria which is the minimum amount of seizures that have to be in the baseline of each patient
        
        14) num_trials:

            (int) -  the number of trials that the endpoints are going to be averaged over

    Outputs:

        1) placebo_RR50_mean:
        
            (float) - the mean of the 50% responder rate of the placebo arm over N (N = num_trials) trials

        2) placebo_RR50_std:
        
            (float) - the standard deviation of the 50% responder rate of the placebo arm over N (N = num_trials) trials
        
        3) placebo_MPC_mean:
        
            (float) - the mean of the median percent change of the placebo arm over N (N = num_trials) trials
        
        4) placebo_MPC_std:

            (float) - the standard deviation of the median percent change of the placebo arm over N (N = num_trials) trials
        
        5) drug_RR50_mean:
        
            (float) - the mean of the 50% responder rate of the drug arm over N (N = num_trials) trials
        
        6) drug_RR50_std:
        
            (float) - the standard deviation of the 50% responder rate of the drug arm over N (N = num_trials) trials
        
        7) drug_MPC_mean:
        
            (float) -  the mean of the median percent change of the drug arm over N (N = num_trials) trials
        
        8) drug_MPC_std:

            (float) - the standard deviation of the median percent change of the drug arm over N (N = num_trials) trials

    '''
    
    # calculate the upper bound on the 95% confidence interval for the normally distributed drug efficacy
    upper_confidence_bound_drug = drug_effect_mu + 1.96*drug_effect_sigma

    # calculate the upper bound on the 95% confidence interval for the normally distributed placebo effect
    upper_confidence_bound_placebo = placebo_effect_mu + 1.96*placebo_effect_sigma
    
    # if the upper bound of the 95% confidence interval for the drug efficacy is less than 1:
    if( upper_confidence_bound_drug <= 1 ):

        # if the upper bound of the 95% confidence interval for the placebo effect is less than 1:
        if( upper_confidence_bound_placebo <= 1):
    
            # calculate the endpoint statistics for the placebo arm
            [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std] = \
                calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                              placebo_effect_mu, placebo_effect_sigma,
                                              drug_effect_mu, drug_effect_sigma, False,
                                              num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                              time_scale_conversion, min_required_baseline_seizure_count, num_trials)

            # calculate the endpoint statistics for the drug arm
            [drug_RR50_mean, drug_RR50_std, drug_MPC_mean, drug_MPC_std] = \
                calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                              placebo_effect_mu, placebo_effect_sigma,
                                              drug_effect_mu, drug_effect_sigma, True,
                                              num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                              time_scale_conversion, min_required_baseline_seizure_count, num_trials)
    
            return [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                    drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std]
        
        # if the upper bound of the 95% confidence interval for the placebo effect is greater than 1:
        else:

            # tell the user that their placebo effect distribution is too likely to generate placebo effects greater than 100%
            raise ValueError('The placebo effect distribution is too likely to generate placebo effects greater than 100% \n(placebo_effect_mu + 1.96*placebo_effect_sigma > 1)')
    
    # if the upper bound of the 95% confidence interval for the drug efficacy is greater than 1:
    else:
        
        # tell the user that their drug efficacy distribution is too likely to generate placebo effects greater than 100%
        raise ValueError('The placebo effect distribution is too likely to generate placebo effects greater than 100% \n(drug_effect_mu + 1.96*drug_effect_sigma > 1)')


# group level parameters of NV model
shape= 24.143
scale = 297.366
alpha = 284.024
beta = 369.628

# parameters determining distribution of placebo effects and drug efficacies
placebo_effect_mu = 0
placebo_effect_sigma = 0
drug_effect_mu = 0.2
drug_effect_sigma = 0.05

# tell the algorithm to generate weekly seizure counts
time_scale_conversion = 7

# RCT design parameters
num_patients_per_arm = 153
num_baseline_intervals = 8
num_testing_intervals = 12

# eligibility criteria
min_required_baseline_seizure_count = 4

# number of trials to calculate endpoint statistics over
num_trials = 5000

# parameters for generating histogram and log-log plot
min_num_weeks = 24
max_num_weeks = 120
num_patients = 10000

start_time_in_seconds = time.time()

# calculate the endpoint statistics under given RCT design parameters
[placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
 drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std    ] = \
    calculate_placebo_and_drug_arm_endpoint_statistics(shape, scale, alpha, beta, 
                                                       placebo_effect_mu, placebo_effect_sigma,
                                                       drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       time_scale_conversion, min_required_baseline_seizure_count, num_trials)

# calculate statistical features of synthetic patient population generated by NV model
[median_monthly_seizure_frequency, log_log_slope, log_log_intercept, r_value,
 monthly_seizure_frequencies, biweekly_log10means, biweekly_log10stds] = \
    get_patient_population_statistical_features(shape, scale, alpha, beta, 
                                                min_num_weeks, max_num_weeks, 
                                                num_patients)


stop_time_in_seconds = time.time()
total_time_in_seconds = stop_time_in_seconds - start_time_in_seconds
total_time_in_minutes = total_time_in_seconds/60
print( 'cpu time: ' +  np.round( str( total_time_in_minutes, 3 ) ) + ' minutes' )

'''
with open( os.getcwd() + '/data.txt' ) as text_file:

    endpoint_responses = np.array([placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                                    drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std    ])

    patient_population_summary_statistics = np.array([median_monthly_seizure_frequency, log_log_slope, 
                                                        log_log_intercept, r_value])

    json.dump(endpoint_responses.tolist(), text_file)
    json.dump(patient_population_summary_statistics.tolist(), text_file)
    json.dump(monthly_seizure_frequencies.tolist(), text_file)
    json.dump(biweekly_log10means.tolist(), text_file)
    json.dump(biweekly_log10stds.tolist(), text_file)
'''


import numpy as np
import scipy.stats as stats
import os
import json
import time


def generate_patient(shape, scale, alpha, beta,
                     num_days_per_interval, num_baseline_intervals, num_total_intervals, min_required_baseline_seizure_count):
    
    acceptable_baseline_rate = False
    
    patient_seizure_counts = np.zeros(num_total_intervals)
    
    while(not acceptable_baseline_rate):
    
        n = np.random.gamma(shape, 1/scale)
        p = np.random.beta(alpha, beta)
    
        mean = n*( (1 - p)/p )
        overdispersion = 1/n
    
        for interval_index in range(num_total_intervals):
        
            rate = np.random.gamma(num_days_per_interval/overdispersion, mean*overdispersion)
            count = np.random.poisson(rate)
        
            patient_seizure_counts[interval_index] = count
        
        if( (num_baseline_intervals != 0) and ( np.sum(patient_seizure_counts[0:num_baseline_intervals]) >= min_required_baseline_seizure_count ) ):
            
            acceptable_baseline_rate = True
        
        elif( num_baseline_intervals == 0 ):
        
            acceptable_baseline_rate = True
    
    return patient_seizure_counts


def get_patient_statistical_features(shape, scale, alpha, beta,
                                                min_num_weeks, max_num_weeks):

    num_testing_intervals = np.random.randint(min_num_weeks, max_num_weeks)

    patient_weekly_seizure_counts = generate_patient(shape, scale, alpha, beta,
                                                      7, 0, num_testing_intervals, 0)
    
    twoWeeks = []
    for y in range(0, len(patient_weekly_seizure_counts)-1, 2):
        twoWeeks.append(patient_weekly_seizure_counts[y]+patient_weekly_seizure_counts[y+1])
    patient_biweekly_seizure_counts = np.array(twoWeeks)

    months = []
    for y in range(0, len(patient_weekly_seizure_counts)-3, 4):
        months.append(patient_weekly_seizure_counts[y]+patient_weekly_seizure_counts[y+1]+patient_weekly_seizure_counts[y+2]+patient_weekly_seizure_counts[y+3])
    patient_monthly_seizure_counts = np.array(months)

    biweekly_log10mean = np.log10(np.mean(patient_biweekly_seizure_counts))
    biweekly_log10std = np.log10(np.std(patient_biweekly_seizure_counts))

    monthly_seizure_frequency = np.mean(patient_monthly_seizure_counts)
    
    return [biweekly_log10mean, biweekly_log10std, monthly_seizure_frequency]


def get_patient_population_statistical_features(shape, scale, alpha, beta, 
                                                min_num_weeks, max_num_weeks, 
                                                num_patients):

    biweekly_log10means = np.zeros(num_patients)
    biweekly_log10stds = np.zeros(num_patients)
    monthly_seizure_frequencies = np.zeros(num_patients)

    for patient_index in range(num_patients):
    
        [biweekly_log10mean, biweekly_log10std, monthly_seizure_frequency] = \
                        get_patient_statistical_features(shape, scale, alpha, beta,
                                                         min_num_weeks, max_num_weeks)
    
        biweekly_log10means[patient_index] = biweekly_log10mean
        biweekly_log10stds[patient_index] = biweekly_log10std
        monthly_seizure_frequencies[patient_index] = monthly_seizure_frequency

    median_monthly_seizure_frequency = np.median(monthly_seizure_frequencies)
    [log_log_slope, log_log_intercept, r_value, _, _] = stats.linregress(biweekly_log10means, biweekly_log10stds)
    
    return [median_monthly_seizure_frequency, log_log_slope, log_log_intercept, r_value, monthly_seizure_frequencies, biweekly_log10means, biweekly_log10stds]


def apply_drug(drug_effect_mu, drug_effect_sigma, patient_seizure_counts,
               num_baseline_intervals, num_testing_intervals):
    
    drug_effect = np.random.normal(drug_effect_mu, drug_effect_sigma)
    
    if(drug_effect >= 1):
        
        drug_effect = 0.99999999
    
    if(drug_effect <= -1):
        
        drug_effect = -0.99999999
        
    testing_counts = patient_seizure_counts[num_baseline_intervals:]
        
    for testing_count_index in range(num_testing_intervals):
        
        num_removed = 0
        testing_count = testing_counts[testing_count_index]
            
        for count_index in range(np.int_(testing_count)):
                
            if(np.random.random() <= np.abs(drug_effect)):
                    
                num_removed = num_removed + np.sign(drug_effect)
            
        testing_counts[testing_count_index] = testing_count  - num_removed
        
    patient_seizure_counts[num_baseline_intervals:] = testing_counts
    
    return patient_seizure_counts


def calculate_endpoints(seizure_counts, num_baseline_intervals, num_patients_per_arm):
    
    baseline_seizure_counts = seizure_counts[:, 0:num_baseline_intervals]
    testing_seizure_counts = seizure_counts[:, num_baseline_intervals:]

    baseline_seizure_frequencies = np.mean(baseline_seizure_counts, 1)
    testing_seizure_frequencies = np.mean(testing_seizure_counts, 1)

    for patient_index in range(num_patients_per_arm):
        if( baseline_seizure_frequencies[patient_index] == 0 ):        
            baseline_seizure_frequencies[patient_index] = 0.00000001

    percent_changes = np.divide(baseline_seizure_frequencies - testing_seizure_frequencies, baseline_seizure_frequencies)
    
    RR50 = 100*np.sum(percent_changes >= 0.5)/num_patients_per_arm
    MPC = 100*np.median(percent_changes)
    
    return [RR50, MPC]


def calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                  drug_effect_mu, drug_effect_sigma, 
                                  num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                  num_days_per_interval, min_required_baseline_seizure_count, num_trials, 
                                  decimal_round):

    num_total_intervals = num_baseline_intervals + num_testing_intervals
    RR50_array = np.zeros(num_trials)
    MPC_array = np.zeros(num_trials)

    for trial_index in range(num_trials):
        
        seizure_counts = np.zeros((num_patients_per_arm, num_total_intervals))

        for patient_index in range(num_patients_per_arm):
        
            seizure_counts[patient_index, :] = generate_patient(shape, scale, alpha, beta,
                                                                num_days_per_interval, num_baseline_intervals, 
                                                                num_total_intervals, min_required_baseline_seizure_count)
        
            seizure_counts[patient_index, :] = apply_drug(drug_effect_mu, drug_effect_sigma, seizure_counts[patient_index, :],
                                                           num_baseline_intervals, num_testing_intervals)

        [RR50, MPC] = calculate_endpoints(seizure_counts, num_baseline_intervals, num_patients_per_arm)
    
        RR50_array[trial_index] = RR50
        MPC_array[trial_index] = MPC

    RR50_mean = np.mean(RR50_array)
    RR50_std = np.std(RR50_array)
    MPC_mean = np.mean(MPC_array)
    MPC_std = np.std(MPC_array)

    return [RR50_mean, RR50_std, MPC_mean, MPC_std]


def calculate_placebo_and_drug_arm_endpoint_statistics(shape, scale, alpha, beta, 
                                                       placebo_effect_mu, drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       num_days_per_interval, min_required_baseline_seizure_count, num_trials, 
                                                       decimal_round):
    
    upper_confidence_bound = drug_effect_mu + 2*drug_effect_sigma
    lower_confidence_bound = drug_effect_mu - 2*drug_effect_sigma
    
    if( (upper_confidence_bound <= 1) and (lower_confidence_bound >= -1) ):
    
        [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std] = \
            calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                          placebo_effect_mu, drug_effect_sigma, 
                                          num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                          num_days_per_interval, min_required_baseline_seizure_count, num_trials, 
                                          decimal_round)

        [drug_RR50_mean, drug_RR50_std, drug_MPC_mean, drug_MPC_std] = \
            calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                          placebo_effect_mu + drug_effect_mu, drug_effect_sigma, 
                                          num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                          num_days_per_interval, min_required_baseline_seizure_count, num_trials, 
                                          decimal_round)
    
        return [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std]
    
    else:
        
        raise ValueError('The 95% confidence interval for the drug effect is either too low or too high.')



shape= 24.143
scale = 297.366
alpha = 284.024
beta = 369.628

placebo_effect_mu = 0
drug_effect_mu = 0.2
drug_effect_sigma = 0.05

num_patients_per_arm = 153
num_days_per_interval = 7
num_baseline_intervals = 8
num_testing_intervals = 12
min_required_baseline_seizure_count = 4
num_trials = 5000
decimal_round = 1

min_num_weeks = 24
max_num_weeks = 120
num_patients = 10000

start_time_in_seconds = time.time()

[placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
 drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std    ] = \
    calculate_placebo_and_drug_arm_endpoint_statistics(shape, scale, alpha, beta, 
                                                       placebo_effect_mu, drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       num_days_per_interval, min_required_baseline_seizure_count, num_trials, 
                                                       decimal_round)

[median_monthly_seizure_frequency, log_log_slope, log_log_intercept, r_value,
 monthly_seizure_frequencies, biweekly_log10means, biweekly_log10stds] = \
    get_patient_population_statistical_features(shape, scale, alpha, beta, 
                                                min_num_weeks, max_num_weeks, 
                                                num_patients)


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

stop_time_in_seconds = time.time()
print((stop_time_in_seconds - start_time_in_seconds)/60)


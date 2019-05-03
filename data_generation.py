import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import textwrap
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


def apply_effect(effect_mu, effect_sigma,
                 patient_seizure_counts,
                 num_baseline_intervals, num_testing_intervals):
    
    effect = np.random.normal(effect_mu, effect_sigma)
    
    if(effect > 1):
        
        effect = 1
    
    testing_counts = patient_seizure_counts[num_baseline_intervals:]
        
    for testing_count_index in range(num_testing_intervals):
        
        num_removed = 0
        testing_count = testing_counts[testing_count_index]
            
        for count_index in range(np.int_(testing_count)):
                
            if(np.random.random() <= np.abs(effect)):
                    
                num_removed = num_removed + np.sign(effect)
            
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
                                  placebo_effect_mu, placebo_effect_sigma,
                                  drug_effect_mu, drug_effect_sigma, drug_arm_flag,
                                  num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                  num_days_per_interval, min_required_baseline_seizure_count, num_trials):

    num_total_intervals = num_baseline_intervals + num_testing_intervals
    RR50_array = np.zeros(num_trials)
    MPC_array = np.zeros(num_trials)

    for trial_index in range(num_trials):
        
        seizure_counts = np.zeros((num_patients_per_arm, num_total_intervals))

        for patient_index in range(num_patients_per_arm):
        
            seizure_counts[patient_index, :] = generate_patient(shape, scale, alpha, beta,
                                                                num_days_per_interval, num_baseline_intervals, 
                                                                num_total_intervals, min_required_baseline_seizure_count)
        
            seizure_counts[patient_index, :] = apply_effect(placebo_effect_mu, placebo_effect_sigma,
                                                            seizure_counts[patient_index, :],
                                                            num_baseline_intervals, num_testing_intervals)
            
            if(drug_arm_flag):

                seizure_counts[patient_index, :] = apply_effect(drug_effect_mu, drug_effect_sigma,
                                                                seizure_counts[patient_index, :],
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
                                                       placebo_effect_mu, placebo_effect_sigma,
                                                       drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       num_days_per_interval, min_required_baseline_seizure_count, num_trials):
    '''

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
        
            (float) - the mean of the normally distributed placebo effect
        
        6) placebo_effect_sigma:

            (float) - the standard deviation of the normally distributed placebo effect
        
        7) drug_effect_mu:
        
            (float) - the mean of the normally distributed drug efficacy
        
        8)  drug_effect_sigma:

            (float) - the standard deviation of the normally distributed drug efficacy

        9) num_patients_per_arm:
        
            (int) - the number of patients for both the placebo arm and the drug arm of a trial
        
        10) num_baseline_intervals:
        
            (int) - the number of intervals in the baseline period of each patient
        
        11) num_testing_intervals:

            (int) - the number of intervals in the testing period of each patient

        12) num_days_per_interval:
        
            (int) -  the number of days per interval in each patients seizure diary
        
        13) min_required_baseline_seizure_count:
        
            (int) - the minimum amount of seizures that have to be in the baseline of each patient
        
        14) num_trials:

            (int) -  the number of trials that the endpoints are going to be averaged over

    Outputs:

        1) placebo_RR50_mean:
        
            (float) - 

        2) placebo_RR50_std:
        
            (float) - 
        
        3) placebo_MPC_mean:
        
            (float) - 
        
        4) placebo_MPC_std:

            (float) -
        
        5) drug_RR50_mean:
        
            (float) -     
        
        6) drug_RR50_std:
        
            (float) - 
        
        7) drug_MPC_mean:
        
            (float) - 
        
        8) drug_MPC_std:

            (float) - 

    '''
    
    # calculate the upper bound on the 95% confidence interval
    upper_confidence_bound = drug_effect_mu + 1.96*drug_effect_sigma
    
    # if the upper bound of the 95% confidence interval is less than 1:
    if( upper_confidence_bound <= 1 ):
    
        [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std] = \
            calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                          placebo_effect_mu, placebo_effect_sigma,
                                          drug_effect_mu, drug_effect_sigma, False,
                                          num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                          num_days_per_interval, min_required_baseline_seizure_count, num_trials)

        [drug_RR50_mean, drug_RR50_std, drug_MPC_mean, drug_MPC_std] = \
            calculate_endpoint_statistics(shape, scale, alpha, beta, 
                                          placebo_effect_mu, placebo_effect_sigma,
                                          drug_effect_mu, drug_effect_sigma, True,
                                          num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                          num_days_per_interval, min_required_baseline_seizure_count, num_trials)
    
        return [placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std]
    
    else:
        
        raise ValueError('The 95% confidence interval for the drug effect is too high.')



shape= 24.143
scale = 297.366
alpha = 284.024
beta = 369.628

placebo_effect_mu = 0
placebo_effect_sigma = 0
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
                                                       placebo_effect_mu, placebo_effect_sigma,
                                                       drug_effect_mu, drug_effect_sigma, 
                                                       num_patients_per_arm, num_baseline_intervals, num_testing_intervals, 
                                                       num_days_per_interval, min_required_baseline_seizure_count, num_trials)

[median_monthly_seizure_frequency, log_log_slope, log_log_intercept, r_value,
 monthly_seizure_frequencies, biweekly_log10means, biweekly_log10stds] = \
    get_patient_population_statistical_features(shape, scale, alpha, beta, 
                                                min_num_weeks, max_num_weeks, 
                                                num_patients)

# add points on log-log plot
fig1 = plt.figure(1)
ax1 = fig1.gca()
plt.scatter(biweekly_log10means, biweekly_log10stds, s = 0.5,  color = 'tab:cyan')
plt.xlim([-0.4, 0.6])
plt.ylim([-.4, 1])

# calculate and add line of best fit on log-log plot
ax1 = plt.gca()
x_vals = np.array(ax1.get_xlim())
y_vals = log_log_intercept + log_log_slope * x_vals
plt.plot(x_vals, y_vals, marker=',', color = 'k', linestyle = '-')

# calculate and add ideal line
x_vals = np.array(ax1.get_xlim())
y_vals = log_log_intercept + 0.7 * x_vals
plt.plot(x_vals, y_vals, color = 'r', linestyle = '--')

# format log-log plot
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
long_title = '2-week seizure count of ' + str( num_patients ) + ' patients'
formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
plt.title(formatted_title, fontsize = 14)
plt.legend(['Line of best fit, slope: ' + str( np.round(log_log_slope, 3) ) + ', R^2 = ' + str( np.round(r_value**2, 3) ), 'target line, slope: 0.7','Individual patient'], fontsize = 12)
plt.xlabel('log10(mean)', fontsize = 14)
plt.ylabel('log10(std.dev)', fontsize = 14)
plt.gray()
fig1.savefig(fname = os.getcwd() + '/Romero-fig1', dpi = 600, bbox_inches = 'tight')

# create histogram
fig2 = plt.figure(2)
[data, bins, _] = plt.hist(monthly_seizure_frequencies, bins='auto', density=True, color = 'tab:cyan')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.axvline(x=median_monthly_seizure_frequency, color='k', linestyle='-')
plt.axvline(x=2.7, color='r', linestyle='--')
plt.legend(['median: ' + str( np.round(median_monthly_seizure_frequency, 2) ), 'target median:2.7'], fontsize = 12)
long_title = 'Histogram of monthly seizure frequencies (' + str( num_patients ) + ' simulated patients)'
formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
plt.title(formatted_title, fontsize = 14)
plt.xlabel('Monthly seizure frequency', fontsize = 14)
long_y_label = 'Fraction of simulated patients'
formatted_y_label = '\n'.join(textwrap.wrap(long_y_label, 30))
plt.ylabel(formatted_y_label, fontsize = 14)
plt.gray()
fig2.savefig(fname = os.getcwd() + '/Romero-fig2', dpi = 600, bbox_inches = 'tight')

simulated_placebo_RR50_mean = np.round(placebo_RR50_mean, decimal_round)
simulated_placebo_RR50_std = np.round(placebo_RR50_std, decimal_round)
simulated_placebo_MPC_mean = np.round(placebo_MPC_mean, decimal_round)
simulated_placebo_MPC_std = np.round(placebo_MPC_std, decimal_round)
simulated_drug_RR50_mean = np.round(drug_RR50_mean, decimal_round)
simulated_drug_RR50_std = np.round(drug_RR50_std, decimal_round)
simulated_drug_MPC_mean = np.round(drug_MPC_mean, decimal_round)
simulated_drug_MPC_std = np.round(drug_MPC_std, decimal_round)

historical_placebo_RR50_mean = 21.1
historical_placebo_RR50_std = 9.9
historical_placebo_MPC_mean = 16.7
historical_placebo_MPC_std = 10.3
historical_drug_RR50_mean = 43.2
historical_drug_RR50_std = 13.1
historical_drug_MPC_mean = 40.9
historical_drug_MPC_std = 11.0

y = [1, 2, 3, 4]

historical_data = [historical_drug_MPC_mean, historical_placebo_MPC_mean,
                       historical_drug_RR50_mean, historical_placebo_RR50_mean]
    
simulated_data = [simulated_drug_MPC_mean, simulated_placebo_MPC_mean,
                    simulated_drug_RR50_mean, simulated_placebo_RR50_mean]
    
historical_std_bars = [historical_drug_MPC_std, historical_placebo_MPC_std,
                        historical_drug_RR50_std, historical_placebo_RR50_std]
    
simulated_std_bars = [simulated_drug_MPC_std, simulated_placebo_MPC_std,
                        simulated_drug_RR50_std, simulated_placebo_RR50_std]
    
ylabels = ['drug arm, MPC', 'placebo arm, MPC', 'drug arm, RR50', 'placebo arm, RR50']
    
xlabel = 'Response percentage'
    
title = 'Efficacy endpoints over simulated and historical RCTs'
    
separation = 0.2
height_const = 0.4
heights = height_const*np.ones(4)
    
[fig4, ax] = plt.subplots()
    
simulated_rects = ax.barh(np.array(y) + separation, simulated_data, height = heights, xerr=simulated_std_bars, label = 'simulated')
historical_rects = ax.barh(np.array(y) - separation, historical_data, height = heights, xerr=historical_std_bars, label = 'historical')
    
i = 0
for rect in simulated_rects:
    text_width = simulated_data[i] + simulated_std_bars[i]
    simumlated_data_string = ('{0:.' + str(decimal_round) + 'f}').format(simulated_data[i])
    simumlated_std_string = ('{0:.' + str(decimal_round) + 'f}').format(simulated_std_bars[i])
    plt.text(1.05*text_width, rect.get_y() + 0.25*rect.get_height(), simumlated_data_string + ' $\pm$ ' + simumlated_std_string)
    i += 1
i = 0
for rect in historical_rects:
    text_width = historical_data[i] + historical_std_bars[i]
    historical_data_string = ('{0:.' + str(decimal_round) + 'f}').format(historical_data[i])
    historical_std_string = ('{0:.' + str(decimal_round) + 'f}').format(historical_std_bars[i])
    plt.text(1.05*text_width, rect.get_y() + 0.25*rect.get_height(), str(historical_data_string) + ' $\pm$ ' + historical_std_string)
    i += 1
    
plt.yticks(y, ylabels, rotation='horizontal')
plt.xlim([0, 100])
plt.xlabel(xlabel)
plt.title(title)
plt.legend()
plt.tight_layout()

plt.show()

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

stop_time_in_seconds = time.time()
print((stop_time_in_seconds - start_time_in_seconds)/60)

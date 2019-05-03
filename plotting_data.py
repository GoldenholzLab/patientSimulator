import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
import json
import textwrap

endpoint_statistics_filename = 'endpoint_statistics'
log_log_histogram_numbers_filename = 'log_log_histogram_numbers'
monthly_seizure_frequencies_filename = 'monthly_seizure_frequencies'
biweekly_log10means_filename = 'biweekly_log10means'
biweekly_log10stds_filename = 'biweekly_log10stds'
decimal_round = 1


with open(os.getcwd() + '/' + endpoint_statistics_filename + '.json', 'r') as text_file:

    endpoint_statistics = np.array(json.load(text_file))

placebo_RR50_mean = endpoint_statistics[0]
placebo_RR50_std  = endpoint_statistics[1]
placebo_MPC_mean  = endpoint_statistics[2]
placebo_MPC_std   = endpoint_statistics[3]
drug_RR50_mean    = endpoint_statistics[4]
drug_RR50_std     = endpoint_statistics[5]
drug_MPC_mean     = endpoint_statistics[6]
drug_MPC_std      = endpoint_statistics[7]

with open(os.getcwd() + '/' + log_log_histogram_numbers_filename + '.json', 'r') as text_file:

    log_log_histogram_numbers = np.array(json.load(text_file))

median_monthly_seizure_frequency = log_log_histogram_numbers[0]
log_log_slope                    = log_log_histogram_numbers[1]
log_log_intercept                = log_log_histogram_numbers[2]
r_value                          = log_log_histogram_numbers[3]


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
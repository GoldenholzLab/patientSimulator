import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
import json
import textwrap

'''

This script plots figures 1, 2, and 4 from the paper associated with this repository. It assumes the existence of

all the needed JSON files containing the relevant information in the same folder as this Python script. If those don;t

exist, generate them with the data_generation.py script.

'''

def load_array_from_json_file(json_file_name):
    '''

    This function loads a 1D Numpy array from a JSON file; the JSON file is assumed to be in the

    same folder as this python script.

    Inputs:

        1) json_file_name:

            (string) - the file name of the JSON file which stores a Numpy array, the JSON file 
                       
                       should be stored in the same folder as this python script

    Outputs:

        1) array:

            (1D Numpy array) - the Numpy array to be loaded from the JSON file

    '''

    # locate the JSON file and load into this python script as an object
    with open(os.getcwd() + '/' + json_file_name + '.json') as json_file:

        # read the Python list from the JSON file and convert it into a 1D Numpy array
        # JSON cannot store Numpy arrays, so NUmpy arrays have to be stored as Python lists
        array = np.array(json.load(json_file))

    return array


def load_data(endpoint_statistics_file_name, 
              log_log_histogram_numbers_file_name, 
              monthly_seizure_frequencies_file_name, 
              biweekly_log10means_file_name, biweekly_log10stds_file_name):
    '''

    This function loads the data needed to fully plot the log-log plot of the L-relationship, the

    histogram of monthly seizure frequencies, and horizontal bar chart of the endpoint response

    statistics.

    Inputs:

        1) endpoint_statistics_file_name:

            (string) - the file name of the JSON file containing the means and standard deviations 

                       of the endpoint responses (RR50, MPC) for both the placebo and drug arms
        
        2) log_log_histogram_numbers_file_name:

            (string) - the file name of the JSON file containing the numerical values needed to 
                       
                       describe the log-log plot as well as the histogram
            
        3) monthly_seizure_frequencies_file_name:

            (string) - the actual data that's needed in order to plot the histogram of monthly 

                       seizure frequencies
        
        4) biweekly_log10means_file_name:

            (string) - the x-axis data needed to plot the log-log plot of the L-relationship
        
        5) biweekly_log10stds_file_name:

            (string) - the y-axis data needed to plot the log-log plot of the L-relationship

    Outputs:

        1) placebo_RR50_mean:

            (float) - the mean of the 50% responder rate of the placebo arm over N (N as determined by user of data_generation.py script) trials

        2) placebo_RR50_std:

            (float) - the standard deviation of the 50% responder rate of the placebo arm over N (N as determined by user of data_generation.py script) trials

        3) placebo_MPC_mean:

            (float) - the mean of the median percent change of the placebo arm over N (N as determined by user of data_generation.py script) trials
            
        4) placebo_MPC_std:

            (float) - the standard deviation of the median percent change of the placebo arm over N (N as determined by user of data_generation.py script) trials
            
        5) drug_RR50_mean:

            (float) - the mean of the 50% responder rate of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        6) drug_RR50_std:

            (float) - the standard deviation of the 50% responder rate of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        7) drug_MPC_mean:

            (float) - the mean of the median percent change of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        8) drug_MPC_std:

            (float) - the standard deviation of the median percent change of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        9) median_monthly_seizure_frequency:

            (float) - the median monthly seizure frequency as estimated from one generated synthetic patient population
            
        10) log_log_slope:

            (float) - the slope of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
            
                      seizure count standard deviation from one generated synthetic patient population
            
        11) log_log_intercept:

            (float) - the intercept of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation from one generated synthetic patient population
            
        12) r_value:

            (float) - the correlation of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation from one generated synthetic patient population
            
        13) biweekly_log10means:

            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count means from one generated synthetic 
        
                               patient population which is used to construct the log-log plot
            
        14) biweekly_log10stds:

            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count standard deviations from one generated 
            
                               synthetic patient population which is used to construct the log-log plot
            
        15) monthly_seizure_frequencies:

            (1D Numpy array) - the array of monthly seizure frequencies from one generated synthetic patient population for 
            
                               which the median monthly seizure frequency was estimated
            
        16) num_patients:

            (int) - the number of patients generated for the syntehtic patient population whose data was used to create the 

                    log-log plot and the histograms
            
    '''

    # load the array containing endpoint response statistics
    endpoint_statistics = load_array_from_json_file(endpoint_statistics_file_name)

    # load the ednpoint response statistics into variables with meaningful variable names
    placebo_RR50_mean = endpoint_statistics[0]
    placebo_RR50_std  = endpoint_statistics[1]
    placebo_MPC_mean  = endpoint_statistics[2]
    placebo_MPC_std   = endpoint_statistics[3]
    drug_RR50_mean    = endpoint_statistics[4]
    drug_RR50_std     = endpoint_statistics[5]
    drug_MPC_mean     = endpoint_statistics[6]
    drug_MPC_std      = endpoint_statistics[7]

    # load the data needed to describe the log-log plot and monthly seizure frequency histogram
    log_log_histogram_numbers = load_array_from_json_file(log_log_histogram_numbers_file_name)

    # load the descriptive data into variables with meaningful variable names
    median_monthly_seizure_frequency = log_log_histogram_numbers[0]
    log_log_slope                    = log_log_histogram_numbers[1]
    log_log_intercept                = log_log_histogram_numbers[2]
    r_value                          = log_log_histogram_numbers[3]

    # load the array of base 10 logarithms of bi-weekly seizure means needed for the x-axis of the log-log plot
    biweekly_log10means = load_array_from_json_file(biweekly_log10means_file_name)

    # load the array of base 10 logarithms of bi-weekly seizure standard deviations needed for the y-axis of the log-log plot
    biweekly_log10stds = load_array_from_json_file(biweekly_log10stds_file_name)
    
    # load the array of monthly seizure frequencies for the histogram
    monthly_seizure_frequencies = load_array_from_json_file(monthly_seizure_frequencies_file_name)
    num_patients  = len(monthly_seizure_frequencies)

    return [placebo_RR50_mean,                placebo_RR50_std,   placebo_MPC_mean,            placebo_MPC_std,
            drug_RR50_mean,                   drug_RR50_std,      drug_MPC_mean,               drug_MPC_std,
            median_monthly_seizure_frequency, log_log_slope,      log_log_intercept,           r_value,
            biweekly_log10means,              biweekly_log10stds, monthly_seizure_frequencies, num_patients]


def plot_l_relationship(log_log_slope, log_log_intercept, r_value, biweekly_log10means, biweekly_log10stds):
    '''

    This function saves a PNG file of the log-log plot of the L-relationship in the same folder as 
    
    this python script.

    Inputs:

        1) log_log_slope:

            (float) - the slope of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
            
                      seizure count standard deviation from one generated synthetic patient population
        
        2) log_log_intercept:

            (float) - the intercept of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation from one generated synthetic patient population

        3) r_value:

            (float) - the correlation of the line of best fit on the log-log plot of biweekly seizure count man vs biweekly 
                      
                      seizure count standard deviation from one generated synthetic patient population

        4) biweekly_log10means:

            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count means from one generated synthetic 
        
                               patient population which is used to construct the log-log plot
            
        
        5) biweekly_log10stds:

            (1D Numpy array) - the array of base 10 logarithms of biweekly seizure count standard deviations from one generated 
            
                               synthetic patient population which is used to construct the log-log plot

    Outputs:

        Technically None

    '''
    # add points on log-log plot
    fig = plt.figure(1)
    ax = fig.gca()
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
    fig.savefig(fname = os.getcwd() + '/Romero-fig2', dpi = 600, bbox_inches = 'tight')


def plot_histogram(median_monthly_seizure_frequency, monthly_seizure_frequencies):
    '''

    This function saves a PNG file of the histogram of the monthly seizure frequencies in the same folder as 

    this Python script.

    Inputs:

        1) median_monthly_seizure_frequency:

            (float) - the median monthly seizure frequency as estimated from one generated synthetic patient population
        
        2) monthly_seizure_frequencies:

            (1D Numpy array) - the array of monthly seizure frequencies from one generated synthetic patient population for 
            
                               which the median monthly seizure frequency was estimated

    Outputs:

        Technically None

    '''

    # create histogram
    fig = plt.figure(2)
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
    fig.savefig(fname = os.getcwd() + '/Romero-fig1', dpi = 600, bbox_inches = 'tight')


def plot_endpoint_responses(placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                            drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std,
                            decimal_round):
    '''

    This function saves a PNG file of the horizontal bar chart of the endpoint response statistics 
    
    in the same folder as this Python script.

    Inputs:

        1) placebo_RR50_mean:

            (float) - the mean of the 50% responder rate of the placebo arm over N (N as determined by user of data_generation.py script) trials

        2) placebo_RR50_std:

            (float) - the standard deviation of the 50% responder rate of the placebo arm over N (N as determined by user of data_generation.py script) trials

        3) placebo_MPC_mean:

            (float) - the mean of the median percent change of the placebo arm over N (N as determined by user of data_generation.py script) trials
            
        4) placebo_MPC_std:

            (float) - the standard deviation of the median percent change of the placebo arm over N (N as determined by user of data_generation.py script) trials
            
        5) drug_RR50_mean:

            (float) - the mean of the 50% responder rate of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        6) drug_RR50_std:

            (float) - the standard deviation of the 50% responder rate of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        7) drug_MPC_mean:

            (float) - the mean of the median percent change of the drug arm over N (N as determined by user of data_generation.py script) trials
            
        8) drug_MPC_std:

            (float) - the standard deviation of the median percent change of the drug arm over N (N as determined by user of data_generation.py script) trial
        
    Outputs:

        Technically None

    '''

    # round the endpoint response statistics 
    simulated_placebo_RR50_mean = np.round(placebo_RR50_mean, decimal_round)
    simulated_placebo_RR50_std = np.round(placebo_RR50_std, decimal_round)
    simulated_placebo_MPC_mean = np.round(placebo_MPC_mean, decimal_round)
    simulated_placebo_MPC_std = np.round(placebo_MPC_std, decimal_round)
    simulated_drug_RR50_mean = np.round(drug_RR50_mean, decimal_round)
    simulated_drug_RR50_std = np.round(drug_RR50_std, decimal_round)
    simulated_drug_MPC_mean = np.round(drug_MPC_mean, decimal_round)
    simulated_drug_MPC_std = np.round(drug_MPC_std, decimal_round)

    # these are the endpoint response statistics from the meta-analysis of 23 RCTs
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
    
    [fig, ax] = plt.subplots()
    
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

    fig.savefig(fname = os.getcwd() + '/Romero-fig4', dpi = 600, bbox_inches = 'tight')


if (__name__ == '__main__'):

    # the file names of the data to be plotted
    endpoint_statistics_file_name = 'endpoint_statistics'
    log_log_histogram_numbers_file_name = 'log_log_histogram_numbers'
    monthly_seizure_frequencies_file_name = 'monthly_seizure_frequencies'
    biweekly_log10means_file_name = 'biweekly_log10means'
    biweekly_log10stds_file_name = 'biweekly_log10stds'

    # the decimal place to which the enpoint response statistics should be rounded to
    decimal_round = 1

    # load the data to be plotted
    [placebo_RR50_mean,               placebo_RR50_std,   placebo_MPC_mean,            placebo_MPC_std,
    drug_RR50_mean,                   drug_RR50_std,      drug_MPC_mean,               drug_MPC_std,
    median_monthly_seizure_frequency, log_log_slope,      log_log_intercept,           r_value,
    biweekly_log10means,              biweekly_log10stds, monthly_seizure_frequencies, num_patients]  = \
        load_data(endpoint_statistics_file_name, 
                  log_log_histogram_numbers_file_name, 
                  monthly_seizure_frequencies_file_name, 
                  biweekly_log10means_file_name, biweekly_log10stds_file_name)

    # plot the log-log plot of the L-relationship
    plot_l_relationship(log_log_slope, log_log_intercept, r_value, biweekly_log10means, biweekly_log10stds)

    # plot the histogram of monthly seizure frequencies
    plot_histogram(median_monthly_seizure_frequency, monthly_seizure_frequencies)

    # plot the horizontal bar char of the endpoint response statistics
    plot_endpoint_responses(placebo_RR50_mean, placebo_RR50_std, placebo_MPC_mean, placebo_MPC_std,
                            drug_RR50_mean,    drug_RR50_std,    drug_MPC_mean,    drug_MPC_std,
                            decimal_round)


'''

This is the main script for this repository. It generates synthetic patients in order to create graphs for the

main paper associated with this article.

'''

from TrialArmClass import TrialArm
import RR50Test
import MPCTest
import math
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import textwrap
import pickle


def checkIntervals(seizureCounts):
    '''
    Checks seizure counts, and returns True if both the mean and standard deviation are non-zero, else returns false. 
    
    This should be one patient's seizure count.
    
    Inputs:
        
        1. seizureCounts
        
            (int) {>=0} - 2D list of seizure counts
            
    Outputs:
        
        1. twoWeeks
        
            returns the list of seizures converted to two week intervals if the mean and standard deviation are non-zero, else returns None
        
        2. months
        
            returns the list of seizures converted to monthly intervals if the mean and standard deviation are non-zero, else returns None
    '''
    # creates array of two-week seizure counts
    twoWeeks = []
    for y in range(0, len(seizureCounts)-1, 2):
        twoWeeks.append(seizureCounts[y]+seizureCounts[y+1])
    
    # creates array of monthly seizure counts
    months = []
    for y in range(0, len(seizureCounts)-3, 4):
        months.append(seizureCounts[y]+seizureCounts[y+1]+seizureCounts[y+2]+seizureCounts[y+3])
    
    # checks the mean and standard deviation of two week slicing to ensure that it is non-zero and non-empty. Monthly seizure counts would likewise be zero/non-zero depending on the two week results, so it does not need to be checked
    if np.mean(twoWeeks) != 0 and np.std(twoWeeks) != 0 and twoWeeks != None and months != None:
        return twoWeeks, months
    return None, None


def plotPatients(patients, image_dir):
    '''
    
    This function takes in all the simulated patient seizure diaries for one placebo arm of a clinical trial,
    
    generates a histogram of monthly seizure frequencies as well as a log-log plot of the mean of the
    
    two-week frequency vs the standard deviation of the two_week frequency based on that data, and saves 
    
    the plots to a specified directory
    
    Inputs:
        
        1) patients:
    
            This parameter is a python list of Patient objects
            
        2) arm:
    
            This parameter is a String that should either be the word 'placebo' or the word 'drug'
            
        3) image_dir:
    
            This parameter is a String which is a path to a directory where the .png images of the 
            
            plots are to be stored
    
    Outputs:
        
        None
    
    '''
    
    # initialize all relevant variables
    numPatients = len(patients)
    mean = []
    stdDev = []
    logMean = []
    logStdDev = []
    twoWeekFreq = []
    currentMonthFreq = []
    averageMonthFreq = []
    
    # calculate the average and standard devation of the monthly frequency, and calculate the mean
    # and standard deviation of the two-week frequency as well as their natural logarithms
    for x in range(numPatients):
        
        twoWeekFreq, currentMonthFreq = checkIntervals(patients[x])
        
        mean.append(np.mean(twoWeekFreq))
        stdDev.append(np.std(twoWeekFreq))
        logMean.append(math.log10(np.mean(twoWeekFreq)))
        logStdDev.append(math.log10(np.std(twoWeekFreq)))
        averageMonthFreq.append(np.mean(currentMonthFreq))
    
    # calculate the slope of the line of best fit between the natural logarithm of the mean
    # and the natural logarithm of the standard deviation
    bestSlope, bestIntercept, bestRValue, bestPValue, bestStdError = stats.linregress(logMean, logStdDev)
    dataMedian = np.median(averageMonthFreq)
        
    #prepare line data for display
    lineDisp = 'slope: ' + str(bestSlope) + ', intercept: ' + str(bestIntercept) + ' \nr-value: ' + str(bestRValue) \
                + ', r squared: ' + str(bestRValue ** 2) + '\np-value: ' + str(bestPValue) + ', std error: ' + str(bestStdError)
    
    print('[median, slope]:\n\n' + str([dataMedian, bestSlope]) + '\n\n' + lineDisp)
    
    # add points on log-log plot
    fig1 = plt.figure(1)
    ax1 = fig1.gca()
    plt.scatter(logMean, logStdDev, s = 0.5,  color = 'tab:cyan')
    #plt.xlim([-.4, 0.6])
    plt.xlim([-0.4, 0.6])
    plt.ylim([-.4, 1])
    
    # calculate and add line of best fit on log-log plot
    ax1 = plt.gca()
    x_vals = np.array(ax1.get_xlim())
    y_vals = bestIntercept + bestSlope * x_vals
    plt.plot(x_vals, y_vals, marker=',', color = 'k', linestyle = '-')
    
    # calculate and add ideal line
    x_vals = np.array(ax1.get_xlim())
    y_vals = bestIntercept + 0.7 * x_vals
    plt.plot(x_vals, y_vals, color = 'r', linestyle = '--')
    
    # format log-log plot
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    long_title = '2-week seizure count of ' + str( numPatients ) + ' patients'
    formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
    plt.title(formatted_title, fontsize = 14)
    plt.legend(['Line of best fit, slope: ' + str( np.round(bestSlope, 3) ) + ', R^2 = ' + str( np.round(bestRValue**2, 3) ), 'target line, slope: 0.7','Individual patient'], fontsize = 12)
    plt.xlabel('log10(mean)', fontsize = 14)
    plt.ylabel('log10(std.dev)', fontsize = 14)
    plt.gray()
    fig1.savefig(fname = image_dir + '/Romero-fig1', dpi = 600, bbox_inches = 'tight')
    
    # create histogram
    fig2 = plt.figure(2)
    [data, bins, _] = plt.hist(averageMonthFreq, bins='auto', density=True, color = 'tab:cyan')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.axvline(x=dataMedian, color='k', linestyle='-')
    plt.axvline(x=2.7, color='r', linestyle='--')
    plt.legend(['median: ' + str( np.round(dataMedian, 2) ), 'target median:2.7'], fontsize = 12)
    long_title = 'Histogram of monthly seizure frequencies (' + str( numPatients ) + ' simulated patients)'
    formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
    plt.title(formatted_title, fontsize = 14)
    plt.xlabel('Monthly seizure frequency', fontsize = 14)
    long_y_label = 'Fraction of simulated patients'
    formatted_y_label = '\n'.join(textwrap.wrap(long_y_label, 30))
    plt.ylabel(formatted_y_label, fontsize = 14)
    plt.gray()
    fig2.savefig(fname = image_dir + '/Romero-fig2', dpi = 600, bbox_inches = 'tight')


def plotPatientWeeks(patient, image_dir):
    '''
    
    This function takes in the simulated patient seizure diary of one patient from the placebo
    
    arm of a clinical trial, and plots the weekly seizure counts for all available weeks
    
    Inputs:
        
        1) patients:
    
            This parameter is a numpy array of weekly seizure counts from one patient
            
        2) image_dir:
    
            This parameter is a String which is a path to a directory where the .png image of the 
            
            plot are to be stored
    
    Outputs:
        
        None
    
    '''    
     # create bar graph of patient's weekly seizure counts
    fig3 = plt.figure(3)
    ax3 = fig3.gca()
    numWeeks = len( patient )
    plt.bar( np.arange( numWeeks ), patient, color='k')
    ax3.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Weeks', fontsize = 14)
    plt.ylabel('weekly seizure counts', fontsize = 14)
    long_title = 'One simulated patient\'s weekly seizure counts over ' + str(numWeeks) + ' weeks'
    formatted_title = '\n'.join(textwrap.wrap(long_title, 40))
    plt.title(formatted_title, fontsize = 14)
    fig3.savefig(fname = image_dir + '/Romero-fig3', dpi = 600, bbox_inches = 'tight')


def generateTrial(shape, scale, alpha, beta, _drugEffect, _placeboEffect, _drugEffectSD, _numScreening, _numBase, _numTestMin, _numTestMax, 
                    _numPlaceboPatients, _numDrugPatients, _screeningMinSzs, _screeningIntervalSize, _screeningIntervalMinSzs, _baseTotalMinSzs,
                    _baseIntervalSize, _baseIntervalMinSzs, _baseSzFree, _uniform):
    '''
    
    This function is essentially a utility function which generates the TrialArm objects corresponding to the placebo arm patients
    
    and the drug arm patients. The code for doing this was organized and collected into this function in order to make the main code
    
    easier to read.
    
    Input:
    
        1. shape
        
            (float) - first group level parameter for the NV model
                
        2. scale
        
            (float) - second group level parameter for the NV model
                
        3. alpha
        
            (float) - third group level parameter for the MV model
                    
        4. beta
        
            (float) - fourth gropu level parameter for the NV model
                    
        5. _drugEffect
        
            (optional) (float) {>=-1 & <=+1} - drug effectiveness; input should be decimal - eg 50% imput as .5
                
                efault: 0 - No effect
                
        6. _drugEffectSD
        
            (optional) (float or int) {>=0} - parameters if drugEffect is to have variability following a normal/gaussian distribution
                
                Default: 0
                
        7. _numScreening
        
            (int) {>=0} - number of weeks to generate for screening seizure counts
                
                Default: 0
        8. _numBase
        
            (optional) (int) {>=0} - number of weeks to generate for baseline seizure counts
                
                Default: 4 - 4 weeks of baseline
        
        9. _numTestMin
            
            (optional) (int) {>=0} - minimum number of weeks to generate for test seizure counts
                
                Default: 5 - miniumum of 5 weeks of testing
        
        10. _numTestMax
        
            (optional) (int) {>=0} - maximum number of weeks to generate for test seizure counts
            
                Default: 5 - maximum of 5 weeks of testing
        
        11. _numPlaceboPatients
            
            (optional) (int) {>=0} - number of patients to generate fro the placebo arm
                
                Default: 100 - 100 patients
                
        12. _numDrugPatients
            
            (optional) (int) {>=0} - number of patients to generate for the drug arm
                
                Default: 100 - 100 patients
        
        13. _screeningMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all screening seizure counts for each patient to be included
                
                Default: 0 - No minimum
        
        14. _screeningIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in screening seizure counts - eg. 4 for evaluating
            
            monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        15. _screeningIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in screening seizure counts necessary for each 
            
            patient to be included. 
            
            The time interval is specified in screeningIntervalSize
                    
                Default: 0 - No minimum
        
        16. _baseTotalMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        17. _baseIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in baseline seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        18. _baseIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        19. _baseSzFree
            
            (int) {} - number of weeks of non-seizure-free periods. None indicates that there is no minimum
                
                Default: None - No minimum
                
        20. _uniform
        
            (bool) {} - True or False statement saying whether or not the lengths of time each patient has their test data
            
                        generated over should be fixed for each patient or distributed uniformly over the min and max number
                        
                        of test weeks specified. If 'uniform' is True, then the fixed length of time each patient has for their
                        
                        test period is the minimum number of test weeks specified
                        
                Default: True - the lengths of time for the test period of each patient is fixed
    
    Outputs:
        
        1) placeboArm:
        
            A TrialArm object corresponding to all the patients in the placebo arm
            
        2) drugArm:
        
            A TrialArm object corresponding to all the patients in the drug arms
    
    '''
    
    # generate the TrialArm object corresponding to all the placebo patients
    placeboArm = TrialArm(shape=shape, scale=scale, alpha=alpha, beta=beta, drugEffect = _placeboEffect, drugEffectSD=_drugEffectSD, 
                     numScreening=_numScreening, numBase=_numBase, numTestMin=_numTestMin, numTestMax=_numTestMax, 
                     numPatients=_numPlaceboPatients, screeningMinSzs=_screeningMinSzs, screeningIntervalSize=_screeningIntervalSize, 
                     screeningIntervalMinSzs=_screeningIntervalMinSzs, baseTotalMinSzs=_baseTotalMinSzs, baseIntervalSize=_baseIntervalSize, 
                     baseIntervalMinSzs=_baseIntervalMinSzs, baseSzFree=_baseSzFree, uniform = _uniform)
    
    # generate the TrialArm object corresponding to all the drug patients
    drugArm = TrialArm(shape=shape, scale=scale, alpha=alpha, beta=beta, drugEffect = _drugEffect + _placeboEffect, drugEffectSD=_drugEffectSD, 
                     numScreening=_numScreening, numBase=_numBase, numTestMin=_numTestMin, numTestMax=_numTestMax, 
                     numPatients=_numDrugPatients, screeningMinSzs=_screeningMinSzs, screeningIntervalSize=_screeningIntervalSize, 
                     screeningIntervalMinSzs=_screeningIntervalMinSzs, baseTotalMinSzs=_baseTotalMinSzs, baseIntervalSize=_baseIntervalSize, 
                     baseIntervalMinSzs=_baseIntervalMinSzs, baseSzFree=_baseSzFree, uniform=_uniform)
    
    return [placeboArm, drugArm]


def sampleMeanAndStandardDevation(shape, scale, alpha, beta, drugEffect, placeboEffect, drugEffectSD, numScreening, numBase, numTestMin, numTestMax,
                                        numPlaceboPatients, numDrugPatients, screeningMinSzs, screeningIntervalSize, screeningIntervalMinSzs, 
                                        baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, uniform, num_trials, data_dir,
                                        genAggValueRound):
    '''
    
    The purpose of this function was to calculate the the sample averages and sample standard deviations of the 50% Responder Rate and the Median
    
    Percent Change for the placebo and drug arms over many trials.
    
    Input:
    
        1. shape
        
            (float) - first group level parameter for the NV model
                
        2. scale
        
            (float) - second group level parameter for the NV model
                
        3. alpha
        
            (float) - third group level parameter for the MV model
                    
        4. beta
        
            (float) - fourth gropu level parameter for the NV model
                    
        5. _drugEffect
        
            (optional) (float) {>=-1 & <=+1} - drug effectiveness; input should be decimal - eg 50% imput as .5
                
                efault: 0 - No effect
                
        6. _drugEffectSD
        
            (optional) (float or int) {>=0} - parameters if drugEffect is to have variability following a normal/gaussian distribution
                
                Default: 0
                
        7. _numScreening
        
            (int) {>=0} - number of weeks to generate for screening seizure counts
                
                Default: 0
        8. _numBase
        
            (optional) (int) {>=0} - number of weeks to generate for baseline seizure counts
                
                Default: 4 - 4 weeks of baseline
        
        9. _numTestMin
            
            (optional) (int) {>=0} - minimum number of weeks to generate for test seizure counts
                
                Default: 5 - miniumum of 5 weeks of testing
        
        10. _numTestMax
        
            (optional) (int) {>=0} - maximum number of weeks to generate for test seizure counts
            
                Default: 5 - maximum of 5 weeks of testing
        
        11. _numPlaceboPatients
            
            (optional) (int) {>=0} - number of patients to generate fro the placebo arm
                
                Default: 100 - 100 patients
                
        12. _numDrugPatients
            
            (optional) (int) {>=0} - number of patients to generate for the drug arm
                
                Default: 100 - 100 patients
        
        13. _screeningMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all screening seizure counts for each patient to be included
                
                Default: 0 - No minimum
        
        14. _screeningIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in screening seizure counts - eg. 4 for evaluating
            
            monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        15. _screeningIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in screening seizure counts necessary for each 
            
            patient to be included. 
            
            The time interval is specified in screeningIntervalSize
                    
                Default: 0 - No minimum
        
        16. _baseTotalMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        17. _baseIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in baseline seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        18. _baseIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        19. _baseSzFree
            
            (int) {} - number of weeks of non-seizure-free periods. None indicates that there is no minimum
                
                Default: None - No minimum
                
        20. _uniform
        
            (bool) {} - True or False statement saying whether or not the lengths of time each patient has their test data
            
                        generated over should be fixed for each patient or distributed uniformly over the min and max number
                        
                        of test weeks specified. If 'uniform' is True, then the fixed length of time each patient has for their
                        
                        test period is the minimum number of test weeks specified
                        
                Default: True - the lengths of time for the test period of each patient is fixed
        
        21) data_dir:
    
            A string which represents a Path to a directory where the text file will be stored
            
        22) genAggValueRound:
        
            An integer describing the number of decimal places to round the aggregate statistic of the RR50 and MPC endpoints to.
    
    Outputs:
        
        1) trials:
    
            (list of tuples of numpy arrays) {} - A list of trials, where each trial is represented as a 2-element tuple, with the
            
                                                  first element being the placebo arm for that trial, and the second element being
                                                  
                                                  the corresponding drug arm
    '''
    
    # initialize the endpoint outputs
    placeboRR50Outputs = np.zeros(num_trials)
    drugRR50Outputs = np.zeros(num_trials)
    placeboMPCOutputs = np.zeros(num_trials)
    drugMPCOutputs = np.zeros(num_trials)
    RR50DrugEffect = np.zeros(num_trials)
    MPCDrugEffect = np.zeros(num_trials)
    RR50PValueSignifiance = np.zeros(num_trials)
    MPCPValueSignifiance = np.zeros(num_trials)
    trials = []
    
    # over a given number of trials
    for i in range(num_trials):
        
        # calculate a specific trial
        [placeboArm, drugArm] = generateTrial(shape, scale, alpha, beta, drugEffect, placeboEffect, drugEffectSD, numScreening, numBase, numTestMin, numTestMax,
                                                numPlaceboPatients, numDrugPatients, screeningMinSzs, screeningIntervalSize, screeningIntervalMinSzs, 
                                                baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, uniform)
        
        trials.append((placeboArm, drugArm))
        
        # calculate the placebo response, drug response, and p-values for that trial
        [placeboRR50Output, drugRR50Output, RR50PValue, placeboMPCOutput, drugMPCOutput, MPCPValue] = calculateRR50AndMPC(placeboArm, drugArm)
        
        # store the responses and p-values
        placeboRR50Outputs[i] = placeboRR50Output
        drugRR50Outputs[i] = drugRR50Output
        placeboMPCOutputs[i] = placeboMPCOutput
        drugMPCOutputs[i] = drugMPCOutput
        RR50DrugEffect[i] = drugRR50Output - placeboRR50Output
        MPCDrugEffect[i] = drugMPCOutput - placeboMPCOutput
        
        # store whether or not the trial was a scuccess accodring to p-value of endpoint
        if(RR50PValue < 0.05):
            RR50PValueSignifiance[i] = 1
        else:
            RR50PValueSignifiance[i] = 0
        if(MPCPValue < 0.05):
            MPCPValueSignifiance[i] = 1
        else:
            MPCPValueSignifiance[i] = 0
    
    #calculate the mean and standard deviations of the responses as well as the percent of significant trials according to the endpoints
    simulated_placebo_RR50_mean = np.round(100*np.mean(placeboRR50Outputs), genAggValueRound)
    simulated_placebo_RR50_std = np.round(100*np.std(placeboRR50Outputs), genAggValueRound)
    simulated_placebo_MPC_mean = np.round(100*np.mean(placeboMPCOutputs), genAggValueRound)
    simulated_placebo_MPC_std = np.round(100*np.std(placeboMPCOutputs), genAggValueRound)
    simulated_drug_RR50_mean = np.round(100*np.mean(drugRR50Outputs), genAggValueRound)
    simulated_drug_RR50_std = np.round(100*np.std(drugRR50Outputs), genAggValueRound)
    simulated_drug_MPC_mean = np.round(100*np.mean(drugMPCOutputs), genAggValueRound)
    simulated_drug_MPC_std = np.round(100*np.std(drugMPCOutputs), genAggValueRound)
    RR50_drug_effect_sample_mean = np.round(np.mean(RR50DrugEffect), genAggValueRound)
    RR50_drug_effect_sample_std = np.round(np.std(RR50DrugEffect), genAggValueRound)
    MPC_drugEffect_sample_mean = np.round(np.mean(MPCDrugEffect), genAggValueRound)
    MPC_drugEffect_sample_std = np.round(np.std(MPCDrugEffect), genAggValueRound)
    percent_significant_RR50 = np.round(np.sum(RR50PValueSignifiance)/num_trials, genAggValueRound)
    percent_significant_MPC = np.round(np.sum(MPCPValueSignifiance)/num_trials, genAggValueRound)
    
    historical_placebo_RR50_mean = 21.08
    historical_placebo_RR50_std = 9.920
    historical_placebo_MPC_mean = 16.66
    historical_placebo_MPC_std = 10.26
    historical_drug_RR50_mean = 43.19
    historical_drug_RR50_std = 13.10
    historical_drug_MPC_mean = 40.89
    historical_drug_MPC_std = 10.97
    
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
        simumlated_data_string = ('{0:.' + str(genAggValueRound) + 'f}').format(simulated_data[i])
        simumlated_std_string = ('{0:.' + str(genAggValueRound) + 'f}').format(simulated_std_bars[i])
        plt.text(1.05*text_width, rect.get_y() + 0.25*rect.get_height(), simumlated_data_string + ' $\pm$ ' + simumlated_std_string)
        i += 1
    i = 0
    for rect in historical_rects:
        text_width = historical_data[i] + historical_std_bars[i]
        historical_data_string = ('{0:.' + str(genAggValueRound) + 'f}').format(historical_data[i])
        historical_std_string = ('{0:.' + str(genAggValueRound) + 'f}').format(historical_std_bars[i])
        plt.text(1.05*text_width, rect.get_y() + 0.25*rect.get_height(), str(historical_data_string) + ' $\pm$ ' + historical_std_string)
        i += 1
    
    plt.yticks(y, ylabels, rotation='horizontal')
    plt.xlim([0, 100])
    plt.xlabel(xlabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    
    fig4.savefig(fname = data_dir + '/Romero-fig4', dpi = 600, bbox_inches = 'tight')
    
    f = open(data_dir + '/RCT Text Data.txt', 'w+')
    
    f.write('\n\nRR50 placebo, mean (standard deviation): ' + str( np.round(simulated_placebo_RR50_mean, 3) ) 
            + ' ( ' + str( np.round(simulated_placebo_RR50_std, 3 ) ) + ' )\n\n')
    f.write('MPC placebo, mean (standard deviation): ' + str( np.round(simulated_placebo_MPC_mean, 3) ) + ' ( ' 
          + str( np.round(simulated_placebo_MPC_std, 3) ) + ' )\n\n')
    f.write('RR50 drug, mean (standard deviation): ' + str( np.round(simulated_drug_RR50_mean, 3)  ) + ' ( ' 
          + str( np.round(simulated_drug_RR50_std, 3) ) + ' )\n\n')
    f.write('MPC drug, mean (standard deviation): ' + str( np.round(simulated_drug_MPC_mean, 3) ) 
          + ' ( ' +str( np.round(simulated_drug_MPC_std, 3) ) + ' )\n\n')
    f.write('RR50 Drug Effect, mean (standard deviation): ' + str( np.round(RR50_drug_effect_sample_mean, 3) ) + ' ( '
          + str( np.round(RR50_drug_effect_sample_std, 3) ) + ' )\n\n')
    f.write('MPC Drug Effect, mean (standard deviation): ' + str( np.round(MPC_drugEffect_sample_mean, 3) ) + ' ( '
          + str( np.round(MPC_drugEffect_sample_std, 3) ) + ' )\n\n')
    f.write('percentage of trials significant according to 50% Response Rate: ' + str(percent_significant_RR50) + '\n\n')
    f.write('percentage of trials significant according to Median Percent Change: ' + str(percent_significant_MPC) + '\n\n')
    
    f.close()
    
    return trials


def calculateRR50AndMPC(placeboArm, drugArm):
    '''
    
    This function takes in all the patient seizure diaries from both the placebo and the
    
    drug arm, and calculates the 50% Responder Rate and Median Percent Change for both arms
    
    as well as corresponding p-values (Fisher Exact Test for the 50% Responder Rate, and 
    
    the Wilcoxon Signed Rank Test for the Median Percent Change)
    
    
    Inputs:
        
        1) placeboArm:
        
            A TrialArm object corresponding to all the patients in the placebo arm
            
        2) drugArm:
        
            A TrialArm object corresponding to all the patients in the drug arm
    
    Outputs:
        
        1) placeboRR50Output:
    
            This is the 50% Responder Rate of all patients from the placebo arm
            
        2) drugRR50Output:
    
            This is the 50% Responder Rate of all patients from the drug arm
            
        3) RR50Pvalue:
    
            This is the p-value as generated by the Fisher Exact test corresponding to the 50% Responder Rates 
            
            of all patients from both arms of the trial
            
        4) placeboMPCOutput:
    
            This is the Median Percent Change of all patients from the placebo arm
            
        5) drugMPCOutput:
    
            This is the Median Percent Change of all patients from the drug arm
        
        6) MPCPValue:
    
            This is the p-value as generated by the Wilcoxon Signed Rank Test corresponding to the Median Percent 
            
            Change of all patients from both arms of the trial 
    
    '''
    
    # extract the baseline placebo, test placebo, baseline drug, and test drug data
    placeboBase = placeboArm.basePatientSzs
    placeboTest = placeboArm.testPatientSzs
    drugBase = drugArm.basePatientSzs
    drugTest = drugArm.testPatientSzs
    
    # calculate the 50% Responder Rate for the placebo arm and the drug arm as well as the corresponding p-value
    contingencyTable = RR50Test.makeTable(placeboBase, placeboTest, drugBase, drugTest)
    placeboRR50Output = contingencyTable[0][0]/(contingencyTable[0][0] + contingencyTable[0][1])
    drugRR50Output = contingencyTable[1][0]/(contingencyTable[1][0] + contingencyTable[1][1])
    RR50PValue = RR50Test.fisherTest(contingencyTable)
    
    # calculate the Median Percent Change for the placebo arm and the drug arm as well as the corresponding p-value
    [placeboPercentChange, drugPercentChange] = MPCTest.patientPercentChange(placeboBase, placeboTest, drugBase, drugTest)
    [placeboMPCOutput, drugMPCOutput] = MPCTest.medianPercentChange(placeboPercentChange, drugPercentChange)
    MPCPValue = MPCTest.wilcoxonTest(placeboPercentChange, drugPercentChange)
    
    return [placeboRR50Output, drugRR50Output, RR50PValue, placeboMPCOutput, drugMPCOutput, MPCPValue]


def plotRR50AndMPC(placeboArm, drugArm, genValueRounding, pValueRounding, image_dir, box):
    
    '''
    
    This function graphs the 50% Responder Rate and the Median Percent Change of all patients from the placebo and drug arms
    
    of a given trial.
    
    Inputs:
        
        1) placeboArm:
        
            A TrialArm object corresponding to all the patients in the placebo arm
            
        2) drugArm:
        
            A TrialArm object corresponding to all the patients in the drug arm
        
        3) genValueRound:
        
            An integer describing the number of decimal places to round the RR50 and MPC endpoints to
         
        4) pValueRounding:
        
            This is how many decimal places the p-value should be rounded to. If the p-value happens to be lower than the
            
            given decimal place, then this function will just say the p-value is lower than the specified decimal place, e.g.,
            
            the function is told to round to 3 decimal places, but the p-value is 0.0000569394, then the function will just 
            
            say p < 0.001.
        
        5) image_dir:
    
            A string which represents a Path to a directory where the .png images of the graphs of the 50% Responder Rate
            
            and the Median Percent Change will be saved to
        
        6) box:
            
    Outputs:
        
        None
    
    '''
    
    # caculate the 50% responder rate and Median Percentage Change
    [placeboRR50Output, drugRR50Output, RR50PValue, placeboMPCOutput, drugMPCOutput, MPCPValue] = calculateRR50AndMPC(placeboArm, drugArm)
    
    # round the p-value for the 50% Responder Rate
    RR50PValueRounded = np.round( RR50PValue, pValueRounding )
    
    # initialize the string to be printed
    RR50PValueString = ''   
    
    # if the p-value is less than the specified decimal place
    if(RR50PValue < 10**(-pValueRounding)):
        # then simply state the p-value as an inequalilty
        RR50PValueString = ' < ' + str( 10**(-pValueRounding) )
    else:
        # otherwise, just print the rounded p-value
        RR50PValueString = ' = ' + str( RR50PValueRounded )
    
    # round the p-value for the 50% Responder Rate
    MPCPValueRounded = np.round( MPCPValue, pValueRounding )
    
    # initialize the string to be printed
    MPCPValueString = ''
    
    #if the p-value is less than the specified decimal place
    if( MPCPValue < 10**(-pValueRounding) ):
         # then simply state the p-value as an inequalilty
        MPCPValueString = ' < ' + str( 10**(-pValueRounding) )     
    else:
        # otherwise, just print the rounded p-value for the Median Percent Change
        MPCPValueString = ' = ' + str( MPCPValueRounded )
    
    # Plot the 50% Responder Rate in the usual way that other RCT papers graph it
    [fig, (ax1, ax2)] = plt.subplots(1, 2)
    rects1 = ax1.bar(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], [100*drugRR50Output, 100*placeboRR50Output], color = 'k')
    ax1.set_xticklabels(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], fontsize = 14)
    ax1.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax1.set_ylim([0, 80])
    ax1.set_ylabel('Percentage of 50% Responders', fontsize = 14)
    longTitle = '50% Responder Rate Results of One Trial (p' + RR50PValueString + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 40))
    ax1.set_title(formatted_title, fontsize = 14)
    for rect in rects1:
        height = rect.get_height()
        ax1.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, genValueRounding), ha = 'center', va = 'bottom', fontsize = 14)
        
    # Plot the Median Percent Change in the usual way that other RCT papers graph it
    rects2 = ax2.bar(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], 
                      [100*drugMPCOutput, 100*placeboMPCOutput], color = 'k')
    ax2.set_xticklabels(['drug (' + str( len(drugArm.allSzs) ) + ' patients)', 'placebo (' + str( len(placeboArm.allSzs) ) + ' patients)'], fontsize = 14)
    ax2.set_yticklabels([0, 10, 20, 30, 40, 50, 60], fontsize = 14)
    ax2.set_ylim([0, 80])
    ax2.set_ylabel('Median Percentage Change', fontsize = 14)
    longTitle = 'Median Percentage Change Results of One Trial (p' + MPCPValueString + ')'
    formatted_title = '\n'.join(textwrap.wrap(longTitle, 45))
    ax2.set_title(formatted_title, fontsize = 14)
    for rect in rects2:
        height = rect.get_height()
        ax2.text(rect.get_x() + rect.get_width()/2, 1.05*height, '%g'%np.round(height, genValueRounding), ha = 'center', va = 'bottom', fontsize = 14)
    
    left  = -0.4  # the left side of the subplots of the figure
    right = 1.3    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.2   # the amount of height reserved for white space between subplots
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    fig.savefig(fname = image_dir + '/Romero-fig5', dpi = 600, bbox_inches = box)


# important meta-parameters
drugEffect = 0.2
drugEffectSD = 0.05
placeboEffect = 0

# randomized clinical trial parameters
numScreening = 0
screeningMinSzs = 0
screeningIntervalSize = 1  
screeningIntervalMinSzs = 0

baseTotalMinSzs = 4
baseIntervalSize = 1
baseIntervalMinSzs = 0
baseSzFree = 0

# group level parameters


shape= 24.143
scale = 297.366
alpha = 284.024
beta = 369.628

'''
shape = 111.313
scale = 296.728
alpha = 296.339
beta = 243.719
'''
'''
shape = 2.802
scale = 50.155
alpha = 59.361
beta = 336.12
'''


############

# remember to change the following line of code when this is ready for publishing

############
# directory to save images to 
image_dir = '/Users/juanromero/Documents/Articles of Interest/Randomized Clinical Trials'

numPlaceboPatients = 10000
numDrugPatients = 1
numBase = 8
numTestMin = 24
numTestMax = 121
uniform = False

[placeboArm, _] = generateTrial(shape, scale, alpha, beta, 0, placeboEffect, 0, numScreening, numBase, numTestMin, numTestMax,
                                        numPlaceboPatients, numDrugPatients, screeningMinSzs, screeningIntervalSize, screeningIntervalMinSzs, 
                                        baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, uniform)

plotPatients(placeboArm.allSzs, image_dir)
plotPatientWeeks(placeboArm.allSzs[0], image_dir)


# number of patients for both arms were taken from trial, and then 50/50 split was provided
numPlaceboPatients = 153
numDrugPatients = 153
numBase = 8
numTestMin = 11
numTestMax = 12
uniform = True
pValueRound = 5
genValueRound = 1
genAggValueRound = 2
num_trials = 5000
box = 'tight'

trials = sampleMeanAndStandardDevation(shape, scale, alpha, beta, drugEffect, placeboEffect, drugEffectSD, numScreening, numBase, numTestMin, numTestMax,
                              numPlaceboPatients, numDrugPatients, screeningMinSzs, screeningIntervalSize, screeningIntervalMinSzs, 
                              baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, uniform,num_trials, image_dir, genAggValueRound)

plotRR50AndMPC(trials[0][0], trials[0][1], genValueRound, pValueRound, image_dir, box)

'''
with open(image_dir + '/histogram_and_log_log_slope_plot_data.pkl', 'wb+') as pkl_file:  
    pickle.dump(placeboArm, pkl_file)
'''
 
with open(image_dir + '/RCT_data.pkl', 'wb+') as pkl_file:
    pickle.dump(trials, pkl_file)

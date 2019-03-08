import numpy as np
from PatientClass import Patient
import math
from scipy import stats

'''

This script generates patients for the parameter optimization scripts, and is not used by the testingTrialResults.py

script.

'''

def generateData(numPatients, shape, scale, alpha, beta):
    '''
    When called, takes in n and p values, along with their associated standard deviations (following a normal/gaussion distribution), 
    
    creates weekly seizure counts of non-treated patients, and returns the median and line of best fit, optionally generating a 
    
    log-log plot of the mean and standard deviation of 2 week seizure frequency, a KDE (Kernel Density Estimate) plot of monthly 
    
    seizure frequencies, and a histogram of monthly seizure frequencies. If more than 20% of patients have a mean/standard deviation 
    
    of 0, this will instead return None, None.
    
    Input:
        1. n
        
            (float or int) {} - µ (mean) of n parameter, parameter of negative binomial generation; n-1 = number of successes
            
        2. p
        
            (float or int) {} - µ (mean) of p parameter, parameter of negative binomial generation where p is the probability of success
            
        3. nSD
        
            (float or int) {} - σ (standard deviation) of n parameter, following a normal/gaussian distribution
            
        4. pSD
        
            (float or int) {} - σ (standard deviation) of n parameter, following a normal/gaussian distribution
            

    Output:
        
        1. dataMedian:
            
            (float) - median of sz freq
            
        2. bestSlope
        
            (float) - slope of line of best fit for the log-log plot
    '''
    
    mean = []
    stdDev = []
    logMean = []
    logStdDev = []
    twoWeekFreq = []
    currentMonthFreq = []
    averageMonthFreq = []
    patientsThrown = 0

    for x in range(numPatients):
        
        '''
        # generates a patient with np.random.uniform(24,576) number of weeks
        # Juan: had to add if-else statement here because of problems with numpy in the O2 cluster
        if(nSD != 0):
            pa = Patient(n=round(abs(np.random.normal(n, abs(nSD)))), p=testingp, numBaseIntervals=np.random.uniform(24,576), numTestIntervals=0, drugEffect=0 )
        else:
            pa = Patient( n=round(n), p=testingp, numBaseIntervals=np.random.uniform(24,576), numTestIntervals=0, drugEffect=0 )
        
        #if x%100 == 0: print('\n', pa.allSzs)
        '''
        
        # changed the number of baseline weeks here
        current_n = np.random.gamma(shape, 1/scale)
        current_p = np.random.beta(alpha, beta)
        pa = Patient(n=current_n, p=current_p, numBaseIntervals=12, numTestIntervals=0, drugEffect=0 )

        
        # if the mean/standard deviation of the two week seizure frequency is non-zero, adds the log10 of the two week seizure frequency
        twoWeekFreq, currentMonthFreq = checkIntervals(pa.allWeeklySzs)
        if ( (twoWeekFreq != None) and (currentMonthFreq != None)  ):
            
            mean.append(twoWeekFreq)
            stdDev.append(twoWeekFreq)
            
            # appends log(mean) and log(std. dev.) of two week slicing to the relevant list
            logMean.append(math.log10(np.mean(twoWeekFreq)))
            logStdDev.append(math.log10(np.std(twoWeekFreq)))
            
            # add average of monthly seizure frequencies for current patient to list of average monthly seizure frequencies
            averageMonthFreq.append(np.mean(currentMonthFreq))
            
        # if either the mean or standard deviation of the two week or monthly slicing is zero, throws away the patient
        else: 
            patientsThrown += 1
        
            # if more than 20% of the patients generated have been thrown away, raise an exception
            if patientsThrown >= 0.2*numPatients:
                return None, None
                #raise ValueError("More than 20% of patients have been thrown away! Please choose new parameters.")

    #calculate line of best fit
    bestSlope, bestIntercept, bestRValue, bestPValue, bestStdError = stats.linregress(logMean, logStdDev)
    
    #Juan: I added this if-else statement because of an error that keeps popping up for me
    if(averageMonthFreq == None or averageMonthFreq == []):
        dataMedian = None
    else:
        dataMedian = np.median(averageMonthFreq)
    
    return dataMedian, bestSlope


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



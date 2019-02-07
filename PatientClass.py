import numpy as np
import seizureMods

'''
defines Patient object, which represents one (1) patient
'''

class Patient():
    '''
    
    Patient object class - Currently, the default generation parameters are for 1 week intervals, and all comments will reflect this
    
    Each patient object will have the following class variables:
    
    Input class variables
        
        1. n
        
            (int) {>=1} - integer parameter of negative binomial generation where n-1 = number of successes
        
        2. p
        
            (float) {>= 0 & <=1} - parameter of negative binomial generation where p = probability of success
        
        3. drugEffect
        
            (float) {>=-1 & <=+1} - % drug effectiveness; input should be decimal - eg 50% effectiveness as .5
        
        4.numBaseIntervals
        
            (int) {>=0} - number of weeks to generate for baseline seizure counts, the default n and p values are weekly seizure counts
        
        5. numTestIntervals
        
            (int) {>=0} - number of weeks to generate for test seizure counts, the default n and p values are weekly seizure counts

    Created class variables
        
        6. allSzs
        
            (int) - 2D matrix of integer test period seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts
        
        7. baseSzs
        
            (int) - 2D matrix of integer test period seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts
        
        8. testSzs
        
            (int) - 2D matrix of integer test period seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts
    '''

    def __init__(self, n=2.25, p=.9, numBaseIntervals=4, numTestIntervals=5, drugEffect=0):
        '''
        
        __init__ () creates template drug-treated patient, given the parameters for the negative binomial, number of baseline/test weeks, 
        
        and drug effectiveness.
        
        possible parameters n, p, numBaseIntervals, numTestIntervals, drugEffect
        
        # n (int) {>0} - parameter of negative binomial generation; n-1 = number of successes
        
        # p (float) {>= 0 & <=1} - parameter of negative binomial generation - p = probability of success    
        
        # numBaseIntervals, numTestIntervals (int) {>=0} - number of weeks to generate for base/drug respectively
        
        # drugEffect (float or int) {-1 to +1} - drug effectiveness; input should be decimal - eg 50% input as .5
        '''
        #variable declaration
        
        self.n = n
        self.p = p

        self.numBaseIntervals = abs(int(numBaseIntervals))
        self.numTestIntervals = abs(int(numTestIntervals))

        # if the absolute value of drug effect is too large, sets it to the largest allowed value
        if drugEffect >=1: self.drugEffect = .999999999
        if drugEffect <=-1: self.drugEffect = -0.999999999
        self.drugEffect = drugEffect%1 if drugEffect >= 0 else (-((abs(drugEffect))%1))

        # The code below generates the baseline seizures for baseline and test time intervals
        [self.rawDailySzs, self.rawWeeklySzs] = self._generateDailyConvertWeekly_(self.n, self.p, self.numBaseIntervals+self.numTestIntervals)
        
        # applies drug (or placebo) effect on the raw seizure counts
        self.allWeeklySzs = seizureMods.applyDrug(self.rawWeeklySzs, self.numTestIntervals, self.drugEffect)

        #takes baseline and test period seizure counts and separates them into different matrices
        self.baseSzs = self.allWeeklySzs[:self.numBaseIntervals]
        self.testSzs = self.allWeeklySzs[self.numBaseIntervals:]


    def _generateDailyConvertWeekly_(self, n, p, numWeeks):
        '''
        
        This function generates daily seizure counts according to a gamma-poisson mixture 
        
        and converts them to weekly seizure counts.
        
        Inputs:
            
            1) n:
    
                parameter to gamma-poisson mixture which must be a real number above zero
                
            2) p:
    
                parameter to gamma-poisson mixture which must be a real number between 0 and 1
                
            3) numWeeks:
    
                length of time in weeks that the seizure count data should cover
                
        Outputs:
            
            1) daily_szs:
    
                numpy array of daily seizure counts
            
            2) weekly_szs:
    
                numpy array of weekly seizure counts, summed up from daily seizure counts
        
        '''
        # generate daily seizure counts with a lenght of time equal the number of weeks multiplied by seven
        daily_szs = np.int_(np.zeros(numWeeks*7))
        
        # for the given number of days
        for i in range(numWeeks*7):
            
            # generate daily seizure outcomes according to gamma-poisson mixture
            rate = np.random.gamma(n, (1-p)/p)
            k = np.random.poisson(rate)
            daily_szs[i] = np.int_(k)
        
        # convert daily seizure counts to weekly seizure counts by summing
        weekly_szs = np.sum( daily_szs.reshape(numWeeks, 7) , 1)
        
        return [daily_szs, weekly_szs]

        
            
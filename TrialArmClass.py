import numpy as np
from PatientClass import Patient
import seizureMods

'''
defines TrailArm object, which represents one (1) trial
'''

class TrialArm(Patient):
    '''
    TrialArm extends Patient - Currently, the default generation parameters are for 1 week intervals, and all comments will reflect this
    
    TrialArm generates a list of curated patients based on input class variables
    
    Class variables
    
    Input class variables
    
        1. shape
        
            (float) - first group level parameter for the NV model
                
        2. scale
        
            (float) - second group level parameter for the NV model
                
        3. alpha
        
            (float) - third group level parameter for the MV model
                    
        4. beta
        
            (float) - fourth gropu level parameter for the NV model
                    
        5. drugEffect
        
            (optional) (float) {>=-1 & <=+1} - drug effectiveness; input should be decimal - eg 50% imput as .5
                
                default: 0 - No effect
                
        6. drugEffectSD
        
            (optional) (float or int) {>=0} - parameters if drugEffect is to have variability following a normal/gaussian distribution
                
                Default: 0
                
        7. numScreening
        
            (int) {>=0} - number of weeks to generate for screening seizure counts
                
                Default: 0
        8.numBase
        
            (optional) (int) {>=0} - number of weeks to generate for baseline seizure counts
                
                Default: 4 - 4 weeks of baseline
        
        9. numTestMin
            
            (optional) (int) {>=0} - minimum number of weeks to generate for test seizure counts
                
                Default: 5 - miniumum of 5 weeks of testing
        
        10. numTestMax
        
            (optional) (int) {>=0} - maximum number of weeks to generate for test seizure counts
            
                Default: 5 - maximum of 5 weeks of testing
        
        11. numPatients
            
            (optional) (int) {>=0} - number of patients to generate
                
                Default: 100 - 100 patients
        
        12. screeningMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all screening seizure counts for each patient to be included
                
                Default: 0 - No minimum
        
        13. screeningIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in screening seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        14. screeningIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in screening seizure counts necessary for each patient to be included. 
            
            The time interval is specified in screeningIntervalSize
                    
                Default: 0 - No minimum
        
        15. baseTotalMinSzs
            
            (int) {} - total minimum number of seizures in the sum of all baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        16. baseIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in baseline seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        17. baseIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        18. baseSzFree
            
            (int) {} - number of weeks of non-seizure-free periods. None indicates that there is no minimum
                
                Default: None - No minimum
                
        19. uniform
        
            (bool) {} - True or False statement saying whether or not the lengths of time each patient has their test data
            
                        generated over should be fixed for each patient or distributed uniformly over the min and max number
                        
                        of test weeks specified. If 'uniform' is True, then the fixed length of time each patient has for their
                        
                        test period is the minimum number of test weeks specified
                        
                Default: True - the lengths of time for the test period of each patient is fixed

    Created class variables
    
        20. allSzs
        
            (int) - 2D matrix of integer all period seizure counts of all patients with the rows as patients and the columns as weekly seizure counts
            
        21. allRawSzs
            
            (int) - 2D matrix of integer all raw (before drug/placebo) seizure counts of all patients with the rows as patients and the columns as weekly seizure counts
            
        22. basePatientSzs
        
            (int) - 2D matrix of integer baseline period seizure counts of all patients with the rows as patients and the columns as weekly seizure counts
        
        23. testPatientSzs
        
            (int) - 2D matrix of integer test period seizure counts of all patients with the rows as patients and the columns as weekly seizure counts
    '''


    def __init__(self, shape, scale, alpha, beta, drugEffect=0, drugEffectSD=0, numScreening=0, numBase=4, numTestMin=5, 
                 
                 numTestMax=5, numPatients=200, screeningMinSzs=0, screeningIntervalSize=1, screeningIntervalMinSzs=0, baseTotalMinSzs=1, 
                 
                 baseIntervalSize=1, baseIntervalMinSzs=0, baseSzFree=None, uniform=True):
        '''
       
        __init__() creates a new TrialArm object, with options to change the following default class variables:
        
        n, nSD, p, pSD, drugEffect, numBase, numTest, numPatients, baseTotalMinSzs, baseIntervalMinSzs
        
        which generate the following class variables:
        
        allSzs, allRawSzs, basePatientSzs, testPatientSzs
        
        '''
        #input cleaning and variable declaration
        self.shape = scale
        self.shape = shape
        self.alpha = alpha
        self.beta = beta
        
        # if the absolute value of drug effect is too large, sets it to the largest allowed value
        if drugEffect >=1: self.drugEffect = .999999999
        if drugEffect <=-1: self.drugEffect = -0.999999999
        self.drugEffect = drugEffect%1 if drugEffect >= 0 else (-((abs(drugEffect))%1))
        self.drugEffectSD = abs(drugEffectSD)
        self.numScreening = round(abs(numScreening))
        self.numBase = round(abs(numBase))
        self.numTestMin = round(abs(numTestMin))
        self.numTestMax = round(abs(numTestMax))
        self.numPatients = round(abs(numPatients))
        self.screeningMinSzs = round(abs(screeningMinSzs))
        self.screeningIntervalSize = round(abs(screeningIntervalSize))
        self.screeningIntervalMinSzs = round(abs(screeningIntervalMinSzs))
        self.baseTotalMinSzs = round(baseTotalMinSzs)
        self.baseIntervalSize = round(abs(baseIntervalSize))
        self.baseIntervalMinSzs = round(baseIntervalMinSzs)
        self.baseSzFree = bool(baseSzFree)
        self.uniform = bool(uniform)
        
        # sets up empty lists
        self.allSzs = []
        self.allRawSzs = []
        self.basePatientSzs = []
        self.testPatientSzs = []

        #generate seizure counts
        for x in range(self.numPatients):
            patientSzs, rawSzs = self.__makePatient()
            self.allSzs.append(patientSzs)
            self.allRawSzs.append(rawSzs)
            self.basePatientSzs.append(patientSzs[:self.numBase])
            self.testPatientSzs.append(patientSzs[self.numBase:])
            patientSzs = []


    def __makePatient(self):
        '''
        __makePatient() is a private subroutine in the TrialArm class that takes no input and 
        
        returns a patient seizure count (both raw and with drug/placebo applied) for the trial 
        
        that has the minimum number of seizure counts outlined in the class variables and with the 
        
        patients generated using class variable values set in __init__ 
        '''
        #checks and returns the patient if it fulfills the minimum seizure requirement, else, generate a new patient 
        while True:
            
            current_n = np.random.gamma(self.shape, 1/self.scale)
            current_p = np.random.beta(self.alpha, self.beta)
            
            if(self.uniform):
                
                pa = Patient(current_n, current_p, self.numBase+self.numScreening, self.numTestMin, np.random.normal(self.drugEffect, self.drugEffectSD))
                
            else:
                
                numTest = np.random.randint(self.numTestMin, self.numTestMax)
                pa = Patient(current_n, current_p, self.numBase+self.numScreening, numTest, np.random.normal(self.drugEffect, self.drugEffectSD))
            
            #checks the base seizure list to see if it meets requirements
            if ( seizureMods.testCriteria(pa.allWeeklySzs, self.numScreening, self.numBase, self.screeningMinSzs, self.screeningIntervalSize,
                                          
                                          self.screeningIntervalMinSzs, self.baseTotalMinSzs, self.baseIntervalSize, self.baseIntervalMinSzs, 
                                          
                                          self.baseSzFree, current_n, current_p) == True ):
                
                return pa.allWeeklySzs, pa.rawWeeklySzs

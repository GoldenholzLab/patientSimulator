import numpy as np

'''
This script holds utility functions for modifying seizure diaries. 

applyDrug reduces (or adds) seizures to the test period depending on the effectiveness of the drug. 

It has an associated private function __applyDrugInternal. testCriteria checks screening and baseline seizure counts to endure that they comply with certain parameters. 

If also has an associated private function __testCriteraInternal.
'''

def applyDrug(szCounts, numTest, drugEffect):
    '''
    applyDrug() takes in seizure counts, the number of test period time units, and the % effectiveness of the drug/placebo, and 
    
    generates a version of the seizure counts where the drug has been applied. In other words, this goes through each patient, 
    
    then loops through each testing interval, and for each seizure, there is a drugEffect % chance of removing (or adding) a seizure.
    
    More than one seizure cannot be added/removed. The function works for both 1D and 2D lists.

    Input:
        
        1. szCounts
        
            (int) {>=0} - 1D matrix or 2D matrix of seizure counts where the rows are patients and the columns time units. 
            
            This can be discontinuous, but if it is a 2D matrix, the first value has to be a 1D matrix and not a single value.
            
        2. numTest
        
            (int) {>=0} - number of unit time units to generate for test seizure counts
            
        3. drugEffect
        
            (float) {>=-1 & <=+1} - % drug effectiveness; input should be decimal - eg 50% effectiveness as .5

    Output:
        
        1. szCounts
        
            (int) {>=0} - 1D matrix or 2D matrix of seizure counts with the drug/placebo applied to it, 
            
            where the rows are patients and the columns time units.
            
    '''
    # checks szCounts as to whether it is a 2D list or not
    if isinstance(szCounts[0], list) == False:

        # applies the drug and returns the new seizure counts
        return __applyDrugInternal(szCounts, numTest, drugEffect)

    else:
        newSzCount = []

        #loops through patients
        for x in range(len(szCounts)):

            # applies the drug and appends the new seizure counts to newSzCount
            newSzCount.append(__applyDrugInternal(szCounts[x], numTest, drugEffect))

        return newSzCount
        

def __applyDrugInternal(szCounts, numTest, drugEffect):
    '''
    
    __applyDrugInternal() is a private function meant to be used exclusively by applyDrug(). 
    
    Given a SINGLE patient's seizurecount, it applies the drug and returns the modified seizure count.
    
    '''
     # loops through each testing interval, and for each seizure, there is a drugEffect % chance of removing (or adding) a seizure. 
     # More than one seizure cannot be added/removed

    # creates a flag to indicate whether to add or remove seizures
    drugPositive = 1 if drugEffect >= 0 else -1

    #loops through test units
    for x in range(numTest):
        numRemoved = 0

        # goes through seizures in one time interval
        for y in range(szCounts[len(szCounts)-numTest+x]):

            # 'rolls the dice' and counts the number of seizures to add/remove
            if (abs(np.random.random())) <= abs(drugEffect):
                numRemoved += 1*drugPositive

        # adds or removes the number of seizures calculated above
        szCounts[len(szCounts)-numTest+x] -= numRemoved
    
    return szCounts


def testCriteria(patientSzs, numScreening, numBase, screeningMinSzs=0, screeningIntervalSize=1, 
                 screeningIntervalMinSzs=0, baseTotalMinSzs=0, baseIntervalSize=1, baseIntervalMinSzs=0, baseSzFree=None, n=1, p=.619):
    '''
    
    testCriteria() takes in seizure counts and the number of weeks in the base and screening, 
    
    with optional arguments as criteria to apply, and returns whether or not the baseline/screening period fit said criteria
    

    Input:
        
        1. patientSzs
        
            (int) {>=0} - 1D matrix or 2D matrix of seizure counts where the rows are patients and the columns weeks. 
            
            This can be discontinuous, but if it is a 2D matrix, the first value has to be a list and not a single value
        
        2. numScreening
        
            (int) {>=0} - number of weeks to generate for screening seizure counts
        
        3. numBase
        
            (int) {>=0} - number of weeks to generate for baseline seizure counts
            
        4. screeningMinSzs
        
            (int) {} - total minimum number of seizures in the sum of all screening seizure counts for each patient to be included
            
                Default: 0 - No minimum
        
        5. screeningIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in screening seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        6. screeningIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in screening seizure counts necessary for each patient to be included. The time interval is specified in screeningIntervalSize
               
                Default: 0 - No minimum
       
        7. baseTotalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in baseline seizure counts necessary for each patient to be included. The time interval is specified in baseIntervalSize
                
                Default: 0 - No minimum
        
        8. baseIntervalSize
            
            (int) {>=1} - size of one unit time interval to be evaluated in baseline seizure counts - eg. 4 for evaluating monthly seizure counts
                
                Default: 1 - 1 week intervals
        
        9. baseIntervalMinSzs
            
            (int) {} - total minimum number of seizures in one unit time interval in baseline seizure counts necessary for each patient to be included
                
                Default: 0 - No minimum
        
        10. baseSzFree
            
            (int) {} - number of weeks of non-seizure-free periods. None indicates that there is no minimum
                
                Default: None - No minimum
        
        11: n
        
            The n parameter used to generate the patient counts for this specific patient (assumes a 1D matrix for patientSzs)
            
        12: p
        
            The p parameter used to generate the patient counts for this specific patient (assumes a 1D matrix for patientSzs)
    Output:
        
        1. _
        
            (boolean) - returns True if patientSzs satisfies all conditions, else returns False

    '''
    
    # checks patientSzs as to whether it is a 2D list or not
    if isinstance(patientSzs[0], list) == False:
        
        # returns True if all flags are true
        return __testCriteriaInternal(patientSzs, numScreening, numBase, screeningMinSzs, screeningIntervalSize,
                                      screeningIntervalMinSzs, baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, n, p)

    else:

        # calls __testCriteriaInternal on each patient
        isTrue = []
        
        for x in range(len(patientSzs)):
            isTrue.append(__testCriteriaInternal(patientSzs[x], numScreening, numBase, screeningMinSzs, screeningIntervalSize, 
                                                 screeningIntervalMinSzs, baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, n, p))

        # returns True if all values in the isTrue list is True
        return all(isTrue)


def __testCriteriaInternal(patientSzs, numScreening, numBase, screeningMinSzs, screeningIntervalSize, 
                           screeningIntervalMinSzs, baseTotalMinSzs, baseIntervalSize, baseIntervalMinSzs, baseSzFree, n, p):
    '''
    
    __testCriteriaInternal() is a private function meant to be used exclusively by testCriteria(). Given a SINGLE patient's seizure count, 
    
    it runs the actual check on patientSzs and returns True if all conditions are satisfied.
    
    '''
    
    # sets screeningTotalMinSzsFlag to true if the total number of screening seizure counts is at least screeningMinSzs
    screeningTotalMinSzsFlag = all(screeningMinSzs<=x for x in patientSzs[:numScreening])

    # ensures that the screeningIntervalSize fits nicely within numScreening, else throws an error
    if numScreening%screeningIntervalSize == 0:
        #variable declaration
        temp = []

        # sets screeningbaseIntervalMinSzsFlag to true if all seizure counts in the specified interval (eg. months) are at least screeningIntervalMinSzs
        for x in range(0,numScreening-(screeningIntervalSize-1),screeningIntervalSize):
            temp.append(sum(patientSzs[x:x+screeningIntervalSize]))
        screeningIntervalMinSzsFlag = all(x>=screeningIntervalMinSzs for x in temp)
    else:
        raise ValueError('The total number of screening time intervals are not divisible by the interval size for criteria - eg. screening size of 5 weeks and expecting seizures per month')
    
    # sets baseTotalMinSzsFlag to true if the total number of baseline seizure counts is at least baseTotalMinSzs
    baseTotalMinSzsFlag = np.sum(patientSzs[numScreening:(numBase+numScreening)])>=baseTotalMinSzs

    # sets baseIntervalMinSzsFlag to true if each baseline interval is at least baseIntervalMinSzs
    baseIntervalMinSzsFlag = all(x>=baseIntervalMinSzs for x in patientSzs[numScreening:(numBase+numScreening)])
    
    # ensures that the baseIntervalSize fits nicely within numBase, else throws an error
    if numBase%baseIntervalSize == 0:
        #variable declaration
        temp = []

        # sets screeningbaseIntervalMinSzsFlag to true if all seizure counts in the specified interval (eg. months) are at least screeningIntervalMinSzs
        for x in range(0,numBase-(baseIntervalSize-1),baseIntervalSize):
            temp.append(sum(patientSzs[x:x+baseIntervalSize]))
        baseIntervalMinSzsFlag = all(x>=baseIntervalMinSzs for x in temp)
    else:
        raise ValueError('The total number of base time intervals are not divisible by the interval size for criteria - eg. base size of 5 weeks and expecting seizures per month')

    

    # sets baseSzFreeFlag to True if the criteria is not applied
    if baseSzFree == None:
        baseSzFreeFlag = True

    # throws an error if the number of seizure-free weeks is larger than possible
    elif baseSzFree > numBase:#numScreening:
        raise ValueError('number of required seizure-free weeks in baseline is greater than the total number of weeks in baseline')
    
    # looks through each section sized baseSzFree and adds True to
    else:
        temp = []

        # looks through each section sized baseSzFree and appends True to a temporary array if there is at least 1 seizure, else appends False
        for x in range(len(patientSzs)-(baseSzFree-1)):
            temp.append(True) if sum(patientSzs[x:x+baseSzFree]>=1) else False
        
        # sets baseSzFreeFlag as True if the entire temporary array is true, else sets it as False
        baseSzFreeFlag = all(temp)
    
    # returns True if all flags are satisfied, else returns False
    return all([screeningTotalMinSzsFlag, screeningIntervalMinSzsFlag, baseTotalMinSzsFlag, baseIntervalMinSzsFlag, 
                baseSzFreeFlag])
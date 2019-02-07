from scipy import stats
import numpy as np

'''
This script holds functions for the median percent change given the baseline and test period seizure counts for both placebo and drug treated arms. 

It can calculate the percent change between a drug and placebo-treated arm (patientPercentChange), 

generate the p-value from the percent change (wilcoxonTest), and calculate the actual meduan percent change (medianPercentChange).
'''

def patientPercentChange(placeboBase, placeboTest, drugBase, drugTest):
    '''
    patientPercentChange() 

    Input:
        1. placeboBase:
            (int) {>=0} - 2D matrix of baseline seizure counts of placebo-treated patients with the rows as patients and the columns as weekly seizure counts
        
        2. placeboTest:
            (int) {>=0} - 2D matrix of test period seizure counts of placebo-treated patients with the rows as patients and the columns as weekly seizure counts
        
        3. drugBase:
            (int) {>=0} - 2D matrix of baseline seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts
        
        4. drugTest
            (int) {>=0} - 2D matrix of test period seizure counts of placebo-treated patients with the rows as patients and the columns as weekly seizure counts

    Output:
        1. placeboPercentChange:
            (float) - list of percent change between baseline and test rates for placebo patients

        2. drugPercentChange
            (float) - list of percent change between baseline rates and test rates for drug-treated patients
    '''
    #make list of % change in placebo trial arm
    placeboPercentChange = []
    for x in range(len((placeboBase))):
        baseRate = sum(placeboBase[x])/len(placeboBase[x])
        testRate = sum(placeboTest[x])/len(placeboTest[x])
        percentChange = (baseRate-testRate)/baseRate if baseRate != 0 else (.0000001-testRate)/.0000001
        placeboPercentChange.append(percentChange)
    
    #make list of % change in drug trial arm
    drugPercentChange = []  
    for x in range(len((drugBase))):
        baseRate = (sum(drugBase[x])/len(drugBase[x]))
        testRate = (sum(drugTest[x])/len(drugTest[x]))
        percentChange = (baseRate-testRate)/baseRate if baseRate !=0 else (.0000001-testRate)/.0000001
        drugPercentChange.append(percentChange)
    
    #return percent changes in placebo and drug arms
    return placeboPercentChange, drugPercentChange


def wilcoxonTest(placeboPercentChange, drugPercentChange):
    '''
    wilcoxonTest() uses the Wilcoxon signed rank test to calculate the median percent change for placebo and drug trial arms, given the percent change lists of percent change from placebo and drug arms
    Input:
        1. placeboPercentChange:
            (float) - list of percent change between baseline and test rates for placebo patients
        2. drugPercentChange
            (float) list of percent change between baseline and test rates for drug-treated patients

    Output:
        1. pvalue:
            (float) - pvalue of % change between placebo and drug
    '''
    #use the Wilcoxon SIgned Rank Test to calculate and return the t value and p value between placebo and drug arms
    _, pvalueOutput = stats.ranksums(placeboPercentChange, drugPercentChange)
    return pvalueOutput


def medianPercentChange(placeboPercentChange, drugPercentChange):
    '''
    medianPercentChange() takes in the patient percent changes and returns the median percent change for both arms
    Input:
        1. placeboPercentChange:
            (float) - list of percent change between baseline and test rates for placebo patients
        2. drugPercentChange
            (float) list of percent change between baseline and test rates for drug-treated patients
    Output:
        1. placeboMPC:
            (float) - median percent change of all placebo patients
        2. drugMPC:
            (float) - median percent change of all drug-treated patients
    '''

    #find the median percent change of placebo and drug
    placeboMPCOutput = np.median(placeboPercentChange)
    drugMPCOutput = np.median(drugPercentChange)
    #return median percent change of placebo and drug
    return placeboMPCOutput, drugMPCOutput
from scipy import stats

'''
This script holds functions for the 50% responder rate given the baseline and test period seizure counts for both placebo and drug treated arms. 

It can make a contingency table (makeTable) of 50% responders, and calculate the p-value using the fisher test (fisherTest).
'''

def makeTable(placeboBase, placeboTest, drugBase, drugTest):
    '''
    makeTable() returns a 2x2 matrix of the 50% response counts for placebo and drug trial arms given matrices of seizure counts of patients in placebo and drug arms

    Input:
        1. placeboBase:
            (int) - 2D matrix of integer baseline seizure counts of placebo-treated patients with the rows as patients and the columns as weekly seizure counts
        
        2. placeboTest:
            (int) - 2D matrix of integer test period seizure counts of placebo-treated patients with the rows as patients and the columns as weekly seizure counts
        
        3. drugBase:
            (int) - 2D matrix of integer baseline seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts
        
        4. drugTest
            (int) - 2D matrix of integer test period seizure counts of drug-treated patients with the rows as patients and the columns as weekly seizure counts

    Output:
        1. contingencyTable
            (2x2 int) - matrix used in Fisher Exact Test, with rows as placebo/drug and columns as a positive/negative responder count respectively
    '''
    #get the number of responders from the placebo base/test matrices
    placeboResponderCount = 0
    for x in range(len((placeboBase))):
        baseRate = (sum(placeboBase[x])/len(placeboBase[x]))
        testRate = (sum(placeboTest[x])/len(placeboTest[x]))
        percentChange = (baseRate-testRate)/baseRate if baseRate !=0 else (.0000001-testRate)/.0000001
        if percentChange>=.5: placeboResponderCount+=1
    
    #get the number of responders from the drug base/test matrices
    drugResponderCount = 0
    for x in range(len((drugBase))):
        baseRate = (sum(drugBase[x])/len(drugBase[x]))
        testRate = (sum(drugTest[x])/len(drugTest[x]))
        percentChange = (baseRate-testRate)/baseRate if baseRate !=0 else (.0000001-testRate)/.0000001
        if percentChange>=.5: drugResponderCount+=1
    
    #creates a contingency table using the number of responders calculated above
    contingencyTable = [[placeboResponderCount, len(placeboBase)-placeboResponderCount], [drugResponderCount, len(drugBase)-drugResponderCount]]
    return contingencyTable
    

def fisherTest(contingencyTable):
    '''
    fisherTest() returns the p value for the 50% responder rate given 

    Input:
        1. contingencyTable (2x2 int list):
            matrix used in Fisher Exact Test, with rows as placebo/drug and columns as a positive/negative responder count respectively

    Output:
        1. pvalue:
            pvalue of placebo vs drug
    '''
    #use the Fisher Exact Test to calculate and return the p value between placebo and drug arms
    _, pvalue = stats.fisher_exact(contingencyTable)
    return pvalue

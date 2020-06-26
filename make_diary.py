import numpy as np


def make_diary(model1or2, number_of_samples, days_per_sample):
    # Use: build one dairy
    # INPUTS:
    #  model1or2 = 1 or 2
    #  number_of_samples = an integer of how many diary samples to produce
    #  days_per_sample = a float of how many days are in each sample. I.e. 7 if weekly, 1 if daily, 0.5 if every 12 hours etc.
    # OUTPUTS:
    #  diary = numpy array of 1 dimension with integer number of seizures
    #  n = the parameter used to generate this patient
    #  p = the other parameter used to generate this patient

    if(model1or2 == 1):

        shape = 24.143
        scale = 297.366
        alpha = 284.024
        beta  = 369.628
    
    elif(model1or2 == 2):

        shape = 111.313
        scale = 296.728
        alpha = 296.339
        beta  = 243.719
    
    else:

        raise ValueError('the \'model1or2\' parameter in the make_diary() function is supposed to take one of two values: \'1\' or \'2\'')

    n = np.random.gamma(shape, 1/scale)
    p = np.random.beta(alpha, beta)

    mean = n*(1 - p)/p
    overdispersion = 1/n

    seizure_diary = np.random.poisson(np.random.gamma(days_per_sample/overdispersion, mean*overdispersion, number_of_samples))

    return [seizure_diary, n, p]


if(__name__=='__main__'):

    model1or2 = 1
    number_of_samples = 3
    days_per_sample = 7

    print(make_diary(model1or2, number_of_samples, days_per_sample))

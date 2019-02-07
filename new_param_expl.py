#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 12:46:55 2019

@author: juanromero
"""

import numpy as np
import pandas as pd
import generation
import time
import argparse


def print4DArray(array_4):
    '''
    
    print4DArray() uses the pandas api to return a string representation of a hierarchical index representation 
    
    of a 4D array. 
    
    Inputs:
        
        1) array_4:
    
            The 4D array to be printed out.
            
    Outputs:
        
       1) array_4_string
    
            A string representation of the 4D array to be printed out
    
    '''
    [dim1, dim2, dim3, dim4] = array_4.shape
    
    array_3_table_list = []
    array_3_keys = dict()
    
    for j in range(dim4):
        
        array_2_table_list  = []
        array_2_keys = dict()
        
        for i in range(dim3):
            
            array_2_table_list.append( pd.DataFrame( array_4[:, :, i, j] ) )
            array_2_keys[i] = array_2_table_list[i]
        
        array_table_3 = pd.concat(array_2_table_list, keys = array_2_keys)
        array_3_table_list.append( array_table_3 )
        array_3_keys[j] = array_3_table_list
        
    array_4_table_list = pd.concat(array_3_table_list, keys = array_3_keys)
    
    array_4_string = array_4_table_list.to_string()
    
    return array_4_string

def lookAroundLocalSpace(central_point, central_point_deltas, decimal_precision):
    '''
    
    lookAroundLocalSpace() takes a 4D point in the parameter space and and generates a local space immediately surrounding
    
    the inputted point which is now the central point of this local space. The space is generated via another 4D point called
    
    the central point step, which dictates how far away or how close the algorithm should look in each dimension in order to 
    
    see the other surrounding points. The function then returns this new local parameter space map.
    
    Inputs:
        
        1) central_point:
    
            A list or tuple of 4 numbers which represents the 4D point (shape, scale, alpha, beta) around which the local space 
            
            will be generated
            
        2) central_point_delta:
    
            A list or tuple of 4 numbers which represents the 4D point (shape_delta, scale_delta, alpha_delta, beta_delta) 
            
            dictating how far or how close the local space will extend in both directions for each dimension
        
        3) decimal_precision:
    
            The number of decimal digits that all the parameter values should be rounded to for each point in each iteration
            
    Outputs:
        
        1) local_space_map:
    
            A 4D numpy array of size (3, 3, 3, 3) which represents the local space surrounding the central point inputted to the
            
            function
        
        2) printable_local_space:
    
            A comprehensible string representing the local space which can be printed to the console
    
    '''
    
    # extract the coordinates of the central point in the parameters space
    shape = central_point[0]
    scale = central_point[1]
    alpha = central_point[2]
    beta = central_point[3]
    
    # extract the magnitude of the step to be taken in each direction from the central point
    shape_delta = central_point_deltas[0]
    scale_delta = central_point_deltas[1]
    alpha_delta = central_point_deltas[2]
    beta_delta = central_point_deltas[3]
    
    # initialize the local parameter space map
    local_space_map = np.array( np.zeros((3, 3, 3, 3)).tolist(), dtype = tuple )
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    local_space_map[i, j, k, l]  = ( np.round( shape + (i - 1)*shape_delta, decimal_precision),
                                                     np.round( scale + (j - 1)*scale_delta, decimal_precision), 
                                                     np.round( alpha + (k - 1)*alpha_delta, decimal_precision),
                                                     np.round( beta + (l - 1)*beta_delta, decimal_precision) ) 
    
    printable_local_space_map = print4DArray(local_space_map)
    
    return [local_space_map, printable_local_space_map]


def probeParameterSpaceLocation(point, num_iter, num_patients_per_iter):
    '''
    
    This function looks at one point in the parameter space and returns the average of the median 
    
    and the average of the slope for a specified number of observations as the estimated median and
    
    estimated slope for that particular point.
    
    Inputs:
        
        1) point:
    
            A tuple of 4 elements, (shape, scale, alpha, beta), which represents a given point in the 
            
            parameter space
            
        2) num_iter:
    
            An integer which is the number of observations needed to calculate the average of the median 
            
            and the average of the slope
        
        3) num_patients_per_iter:
                        
            An integer which is the number of patients generated per observation
            
    Outputs:
        
        1) median_average:
    
            The average of the median for the given location in the parameter space, with the number of
            
            observations used to calculate this average specified by the input parameter num_iter
            
        1) slope_average:
    
            The average of the slope for the given location in the parameter space, with the number of
            
            observations used to calculate this average specified by the input parameter num_iter
    
    '''
    
    # extract the given parametes from the given point
    shape = point[0]
    scale = point[1]
    alpha = point[2]
    beta = point[3]
    
    # initialize the lists of medians and slopes to average over
    median_list = []
    slope_list = []
    
    # initialize flag for whether or not point can be probed
    could_be_probed = True
    
    # initialize counter for how many times an observation couldn't be generated
    observations_not_generated = 0
    
    # initialize counter for while loop
    counter = 0
    
    # while the counter is below the specified number of observations and the point is not yet determined to be intractable
    while((counter < num_iter) and (could_be_probed)):
        
        # try to generate an observation
        [current_median, current_slope] = generation.generateData(num_patients_per_iter, shape, scale, alpha, beta)
            
        # if this attempt failed
        if((current_median == None) or (current_slope == None)):
            
            # increase the number of failed observations by one
            observations_not_generated += 1
            
            # check to see if the number of failed observations is equal to the specified number of observations
            if(observations_not_generated >= num_iter):
                
                # if so, say that the point is intractable
                could_be_probed = False
                
        # otherwise
        else:        
            # append the observed median and slope to their respective lists
            median_list.append(current_median)
            slope_list.append(current_slope)
            
            # iterate the counter
            counter += 1
    
    # if there were no problems in probing the given point
    if(could_be_probed):
        
        # take the respective averages of the medians and slopes over all the observations
        median_average = np.mean(median_list)
        slope_average = np.mean(slope_list)
        
        # and return both of those averages
        return [median_average, slope_average]
    
    else:
        
        # otherwise, return not a number for both values
        return [np.nan, np.nan]


def probeParameterSpaceLocationError(point, num_iter, num_patients_per_iter, target_values):
    '''
    
    This function calculates the error for any given point using the target median, the
    
    target slope, and the estimated median as well as the estimated slope. The error is calculated
    
    via the square of the L-2 norm along the median and slope.
    
    Inputs:
        
        1) point:
    
            A tuple of 4 elements, (shape, scale, alpha, beta), which represents a given point in the 
            
            parameter space
            
        2) num_iter:
    
            An integer which is the number of observations needed to calculate the average of the median 
            
            and the average of the slope
            
        3) num_patients_per_iter:
                        
            An integer which is the number of patients generated per observation
            
        4) target_values:
    
            A 2-element list containing the median and slope that from the desired probability distribution
            
    Outputs:
        
        1) point_error:
    
            The final computed error for the given point =, which is calculated as such below:
                
                
                            point_error = (median_error)^2 + 10*(slope_error)^2
                            
            
            where median_error and slope_error are defined as:
                
                
                                median_error = median_average - target_median
                
                                slope_error = slope_average - target_slope
    '''
    
    target_median = target_values[0]
    target_slope = target_values[1]
    
    # calculate the average median and average slope for the given point
    [median_average, slope_average] = probeParameterSpaceLocation(point, num_iter, num_patients_per_iter)
    
    # calculate the individual error for the median and indivdual error for the slope
    median_error = median_average - target_median
    slope_error = slope_average - target_slope
    
    # calculate the error for the entire point based off of the L-2 norm combining the 
    # error for the median and the error for the slope
    point_error = np.square(median_error) + 100*np.square(slope_error)
    
    return point_error


def checkParameterBoundaries(parameter_point):
    '''
    
    checkParameterBoundaries() looks at shape, scale, alpha, and beta to see if any of them are zero. 
    
    If one or more of those conditions aren't true, then the function returns a False boolean value.
    
    Otherwise, it returns a True boolean value.
    
    Inputs:
        
        1) parameter_point:
                        
            A list or tuple of 4 numbers which represents the 4D point in the parameter space which needs to be checked
    
    Outputs:
        
        1) The output is a simple boolean value, it's either True or False
    
    '''
    
    # extract shape, scale, alpha, and beta from the given parameter point
    shape = parameter_point[0]
    scale = parameter_point[1]
    alpha = parameter_point[2]
    beta = parameter_point[3]
    
    # calculate the boolean values on whether or not the boundaries of the parameters are still true
    bool1 = shape > 0
    bool2 = scale > 0
    bool3 = alpha > 0
    bool4 = beta > 0
    bool5  = bool1 and bool2 and bool3 and bool4
    
    # return True if the boundaries were not violated, return False otherwise
    if(bool5):
        return True
    else:
        return False


def probeLocalParameterSpaceError(local_space_map, num_iter, num_patients_per_iter, target_values, decimal_precision):
    '''
    
    This function estimates the error at each point within the local parameter space, and returns a map
    
    of all the estimated errors at each point within the local space.
    
    Inputs:
        
        1) local_space_map:
    
            A 4D Numpy array representing the local space of all the points immediately surrounding the 
            
            central point; this parameter shoudl be the output of the function lookAroundLocalSpace()
            
        2) num_iter:
    
            An integer which is the number of observations needed to calculate the average of the median 
            
            and the average of the slope
        
        3) num_patients_per_iter:
                        
            An integer which is the number of patients generated per observation
            
        4) target_values:
    
            A 2-element list containing the median and slope that from the desired probability distribution
            
        5) decimal_precision:
    
            An integer which decides how many decimal places to round each number in the output to
            
    Outputs:
        
        1) local_space_errors:
    
            A 4D Numpy array consisting of the errors for each point in the local_space parameter
        
        2) printable_local_space_error:
    
            A comprehensible string representing the local space of errors which can be printed to the console
    
    '''
    
    # initialize the array representing the local space of parameters
    local_space_errors = np.zeros((3, 3, 3, 3))
    
    # for each in the local space
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
            
                    if( checkParameterBoundaries(local_space_map[i, j, k, l]) ):
            
                        # calculate the error for each point in the local space
                        local_space_errors[i, j] = \
                                                    \
                        probeParameterSpaceLocationError(local_space_map[i, j, k, l], num_iter, num_patients_per_iter, target_values)
            
                    else:
            
                        local_space_errors[i, j, k, l] = np.nan
            
    
    # generate a printable string version of the errors in the local parameter space
    printable_local_space_errors = print4DArray( np.round( local_space_errors, decimal_precision ) )
    
    return [local_space_errors, printable_local_space_errors]


def findMinimalErrorPoint(local_space_map, local_space_errors):
    '''
    
    This function finds the point within the local space that has the minimum error as well
    
    as what the error is.
    
    Inputs:
        
        1) local_space_map:
    
            A 4D Numpy array representing the local space of all the points immediately surrounding the 
            
            central point; this parameter should be the output of the function lookAroundLocalSpace()
            
        2) local_space_errors:
    
            A 4D Numpy array consisting of the errors for ecah point in the local_space parameter; 
            
            this parameter should be the output of the function probeLocalParameterSpaceError()
            
    Outputs:
        
        1) min_error_point:
    
            The point in the local parameter space with the lowest error, taking the form:
                
                                    (shape, scale, alpha, beta)
            
        2) min_error:
    
            The error corresponding to min_error_point
    
    '''
    
    # find the indices of the point in the space of local parameters with the lowest error
    min_error_point_index = np.unravel_index( np.nanargmin(local_space_errors, axis=None), local_space_errors.shape )
    
    # find the shape,scale, alpha, and beta coordinates of the the point with the lowest error
    min_error_shape_index = min_error_point_index[0]
    min_error_scale_index = min_error_point_index[1]
    min_error_alpha_index = min_error_point_index[2]
    min_error_beta_index = min_error_point_index[3]
    
    # find the minimal error point and its corresponding error
    min_error_point = local_space_map[min_error_shape_index, min_error_scale_index, min_error_alpha_index, min_error_beta_index]
    min_error = local_space_errors[min_error_shape_index, min_error_scale_index, min_error_alpha_index, min_error_beta_index]
    
    return [min_error_point, min_error]


def step(central_point, central_point_deltas, decimal_precision, num_iter, num_patients_per_iter, 
         target_values, journey_flag):
    '''
    
    This function takes a central point, looks at the local space all around the central point, finds the
    
    next point with the lowest error, and returns that point along with its associated error.
    
    Inputs:
        
        1) central_point:
    
            A tuple of four postivie integers, taking the form:
                
                                        (shape, scale, alpha, beta)
        
        2) central_point_deltas:
    
            A tuple of four positive floating_point numbers, taking the form:
                
                            (shape_delta, scale_delta, alpha_delta, beta_delta)
        
            
        3) decimal_precision:
    
            An integer which decides how many decimal places to round each number in the output to
            
        4) num_iter:
    
            An integer which is the number of observations needed to calculate the average of the median 
            
            and the average of the slope
            
        5) num_patients_per_iter:
                        
            An integer which is the number of patients generated per observation
            
        6) target_values:
    
            A 2-element list containing the median and slope that from the desired probability distribution
        
        7) journey_flag:
    
            A boolean value that's True if the algorithm should output its location at each iteration,
            
            and False otherwise
            
        
        Outputs:
            
            1) new_central_point:
    
                A tuple with the same form as the input parameter central_point, but with different values.
                
                This output represents the new point that the algorithm needs to move to.
                
            2) min_error:
    
                This is a floating-point number which is the error calculated at new_central_point
    
    '''
    
    # look around the local space and find the point with the miinimum error
    [local_space_map, printable_local_map] = lookAroundLocalSpace(central_point, central_point_deltas, decimal_precision)
    [local_space_errors, printable_local_space_errors] = \
        probeLocalParameterSpaceError(local_space_map, num_iter, num_patients_per_iter, target_values, decimal_precision)
    [min_error_point, min_error] = findMinimalErrorPoint(local_space_map, local_space_errors)
    
    # print out the data if journey_flag is True
    if(journey_flag):
        
        print('local space of parameters:\n\n' +  printable_local_map + '\n\n\n\nlocal space of parameter errors:\n\n' 
              
              + printable_local_space_errors + '\n\n\n\n[minimum error point, minimum error]:\n\n' 
              
              + str([min_error_point, min_error]) + '\n\n\n\nmedian estimate, slope estimate]:\n\n' + 
        
              str( probeParameterSpaceLocation( min_error_point, num_iter, num_patients_per_iter) ) + '\n\n')
    
    # set minimal error point as new central point
    new_central_point = (min_error_point[0], min_error_point[1], min_error_point[2], min_error_point[3])
    
    return [new_central_point, min_error]


def init_central_point(init_param_boundaries, decimal_precision):
    '''
    
    init_central_point() randomly chooses an initial point residing in the 4D parameter space within the bounds set by the inputs.
    
    Inputs:
        
      1) init_param_boundaries:
    
            A list of 4 sub-lists, with each sublist being a 2-element of floats that is the lower and upper
            
            boundary for shape, scale, alpha, and beta parameters, respectively
      
      2) decimal_precision:
                        
            The number of decimal digits that all the initial values should be rounded to
            
    Outputs:
        
        1) central_point:
    
            A list or tuple of 4 numbers which represents the intial 4D point (shape, scale, alpha, beta) that will be given to the
            
            rest of the algorithm
    '''
    shape_lower_Boundary_init = init_param_boundaries[0][0]
    shape_upper_Boundary_init = init_param_boundaries[0][1]
    scale_lower_Boundary_init = init_param_boundaries[1][0]
    scale_upper_Boundary_init = init_param_boundaries[1][1]
    alpha_lower_Boundary_init = init_param_boundaries[2][0]
    alpha_upper_Boundary_init = init_param_boundaries[2][1]
    beta_lower_Boundary_init = init_param_boundaries[3][0]
    beta_upper_Boundary_init = init_param_boundaries[3][1]
    
    shape = np.round(np.random.uniform(shape_lower_Boundary_init, shape_upper_Boundary_init), decimal_precision)
    scale = np.round(np.random.uniform(scale_lower_Boundary_init, scale_upper_Boundary_init), decimal_precision)
    alpha = np.round(np.random.uniform(alpha_lower_Boundary_init, alpha_upper_Boundary_init), decimal_precision)
    beta = np.round(np.random.uniform(beta_lower_Boundary_init, beta_upper_Boundary_init), decimal_precision)
    
    central_point = (shape, scale, alpha, beta)
    
    return central_point


def algorithm(init_param_boundaries, central_point_deltas, decimal_precision, num_iter, num_patients_per_iter, 
              max_iter, target_values, max_tolerable_error, journey_flag, output_file_name):
    '''
    
    This function is the function which implements the entire algorithm.
    
    Inputs:
        
        1) init_param_boundaries:
    
            A list of 4 sub-lists, with each sublist being a 2-element of floats that is the lower and upper
            
            boundary for shape, scale, alpha, and beta parameters, respectively
        
        2) central_point_deltas:
    
            A tuple of four positive floating_point numbers, taking the form:
                
                            (shape_delta, scale_delta, alpha_delta, beta_delta)
        
            
        3) decimal_precision:
    
            An integer which decides how many decimal places to round each number in the output to
            
        4) num_iter:
    
            An integer which is the number of observations needed to calculate the average of the median 
            
            and the average of the slope
        
        5) num_patients_per_iter:
                        
            An integer which is the number of patients generated per observation
            
        6) max_iter:
    
            The maximum number of iterations allowed for the algorithm to operate under
            
        6) target_values:
    
            A 2-element list containing the median and slope that from the desired probability distribution
            
        7) max_tolerable_error:
    
            The maximum error that the algorithm must reach before it can stop
        
        8) journey_flag:
    
            A boolean value that's True if the algorithm should output its location at each iteration,
            
            and False otherwise
            
        9) output_file_name:
    
            A string which is the nameof the file to which the output will be written, the output being the set of
            
            parameters that the algorithm stops at
            
    Outputs:
        
        None
    
    '''
    
    # initialize the number of steps taken by the algorithm
    num_steps = 0
    
    # initialize the current minimum error at a ridiculously high number
    min_error = 2000000000
    
    # randoml initialize the central point within certain pre-specified boundaries
    central_point = init_central_point(init_param_boundaries, decimal_precision)
    
    # initialize flag for whether or not a local minima was hit
    local_minima = False
    
    # while the error at the current point is higher than the maximum tolerable error
    while( (min_error > max_tolerable_error) and (not local_minima) and (num_steps < max_iter) ):
        
        if(journey_flag):
            
            print('-----------------------------------------')
            print('step #: ' + str(num_steps) + '\n\n')
        
        start_time_sec = time.time()
        
        # find the neighboring point within the local space that has the current lowest error
        [new_central_point, min_error] = \
                step(central_point, central_point_deltas, decimal_precision, num_iter, num_patients_per_iter, 
                     target_values, journey_flag)
        
        end_time_sec = time.time()
        total_time_sec = end_time_sec  - start_time_sec
        total_time_min = total_time_sec/60
        
        # check to see whether or not the algorithm actually moved or not                        
        if( new_central_point == central_point ):
            
            # if so, then raise flag
            local_minima = True
        
        # move to that point and make it the new central point
        central_point = new_central_point
        
        # iterate the number of steps taken
        num_steps += 1
        
        if(journey_flag):
            
            print('total time (minutes): ' + str(total_time_min) + '\n-----------------------------------------')
    
    if(local_minima):
        print('\n\nlocal minimum reached')
    else:
        print('\n\nreached end of computation')
    
    with open(output_file_name, 'w+') as output_file:
        
        output_file.write( 'minimum cost: ' + str(min_error) + ', point: ' + str(central_point) )



if (__name__ == '__main__'):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('array', nargs = '+')
    args = parser.parse_args()
    arg_array = args.array
    
    shape_lower_Boundary_init = float(arg_array[0])
    shape_upper_Boundary_init = float(arg_array[1])
    scale_lower_Boundary_init = float(arg_array[2])
    scale_upper_Boundary_init = float(arg_array[3])
    alpha_lower_Boundary_init = float(arg_array[4])
    alpha_upper_Boundary_init = float(arg_array[5])
    beta_lower_Boundary_init = float(arg_array[6])
    beta_upper_Boundary_init = float(arg_array[7])
    
    shape_delta = float(arg_array[8])
    scale_delta = float(arg_array[9])
    alpha_delta = float(arg_array[10])
    beta_delta = float(arg_array[11])
    
    target_median = float(arg_array[12])
    target_slope = float(arg_array[13])
    num_iter = int(arg_array[14])
    max_iter = int(arg_array[15])
    num_patients_per_iter = int(arg_array[16])
    decimal_precision = int(arg_array[17])
    
    max_tolerable_error_dB = float(arg_array[18])
    max_tolerable_error = np.power(10, max_tolerable_error_dB/10)
    
    if(arg_array[19] == 'True'):
        
        journey_flag = True
    
    elif(arg_array[19] == 'False'):
        
        journey_flag = False
    
    else:
        
        raise ValueError('20th argument to new_param_expl.py needs to be a boolean value')
        
    output_file_name = arg_array[20]
    
    ''''
    shape_lower_Boundary_init = 10
    shape_upper_Boundary_init = 50
    scale_lower_Boundary_init = 180
    scale_upper_Boundary_init = 220
    alpha_lower_Boundary_init = 140
    alpha_upper_Boundary_init = 180
    beta_lower_Boundary_init = 80
    beta_upper_Boundary_init = 120
    
    shape_delta = 1
    scale_delta = 1
    alpha_delta = 1
    beta_delta = 1
    
    target_median = 2.7
    target_slope = 0.7
    num_iter = 30
    max_iter = 1
    num_patients_per_iter = 1000
    decimal_precision = 3
    
    max_tolerable_error_dB = -30
    max_tolerable_error = np.power(10, max_tolerable_error_dB/10)
    
    journey_flag = True
    output_file_name = 'simulation_output.txt'
    '''
    
    central_point_deltas = (shape_delta, scale_delta, alpha_delta, beta_delta)
    target_values = [target_median, target_slope]
    init_param_boundaries = [[shape_lower_Boundary_init, shape_upper_Boundary_init],
                             [scale_lower_Boundary_init, scale_upper_Boundary_init],
                             [alpha_lower_Boundary_init, alpha_upper_Boundary_init],
                             [beta_lower_Boundary_init, beta_upper_Boundary_init]]
    
    algorithm(init_param_boundaries, central_point_deltas, decimal_precision, num_iter, num_patients_per_iter, 
              max_iter, target_values, max_tolerable_error, journey_flag, output_file_name)


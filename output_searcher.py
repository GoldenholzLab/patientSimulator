#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 09:01:54 2019

@author: juanromero
"""

import glob
import numpy as np


output_directory = '/n/scratch2/jmr95/parameter_optimization'
file_names = '/*.txt'

output_file_names = output_directory + file_names
files = glob.glob(output_file_names)

output_costs = []
output_points = []

for file_name in files:
    try:
        with open(file_name) as file:
            
            output_string = file.read()
            output_strings = output_string.split(", ")
        
            output_cost_strings = output_strings[0].split(': ')
            output_cost = float(output_cost_strings[1])
            output_shape_strings = output_strings[1].split('(')
            output_beta_strings = output_strings[4].split(')')
            
            output_shape = float(output_shape_strings[1])
            output_scale = float(output_strings[2])
            output_alpha = float(output_strings[3])
            output_beta = float(output_beta_strings[0])
            
            output_costs.append(output_cost)
            output_points.append((output_shape, output_scale, output_alpha, output_beta))
            
            '''
            output_point_string = '(' + str(output_shape) + ', ' + str(output_scale) + ', ' + \
                                    str(output_alpha) + ', ' + str(output_beta) + ')'
            
            print('\n\n' + str(output_cost) + '\n' + output_point_string)
            '''
            
    except IOError as exc:
            
        print(exc)
            
min_cost_index = np.nanargmin( np.array(output_costs) )
        
print('\n\n' + str(output_costs[min_cost_index]) + '\n' + str(output_points[min_cost_index]))            

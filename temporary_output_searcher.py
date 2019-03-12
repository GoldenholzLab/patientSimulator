import glob
import numpy as np

output_directory = '/home/jmr95/patientSimulator'
error_file_names = '/*.err'

error_file_names = output_directory + error_file_names
output_file_names = []
error_file_names = glob.glob(error_file_names)
costs = []

for error_file_name in error_file_names:
    try:
        with open(error_file_name, 'r') as error_file:
            if( error_file.readline() == 'srun: Job step aborted: Waiting up to 62 seconds for job step to finish.\n'):
                [file_name, file_extension] = error_file_name.split('.')
                output_file_name = file_name + '.out'
                output_file_names.append(output_file_name)
    except IOError as exc:
        print(exc)

for output_file_name in output_file_names:
    try:
        with open(output_file_name, 'r') as output_file:
            lines = output_file.readlines()
            lines.reverse()
            #print(lines[0])
            #print(lines[1])
            #print(lines[2])
            #print(lines[3])
            print(lines[4])
            #print(lines[5])
            '''
            if(lines[3] == 'median estimate, slope estimate]:'):
                lines = lines[0:10]
                lines.reverse()
                line1 = lines[2]
                #line2 = lines[8]
                #print(output_file_name )
                #print(line1.split(','))
                cost = float(line1.split(',')[4].split(']')[0])
                costs.append(cost)
            '''
    except IOError as exc:
        print(exc)

#print(np.argmin(costs))
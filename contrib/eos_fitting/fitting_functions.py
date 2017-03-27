import os
import sys
import numpy as np

def read_fitting_file(filename):
    with open(filename, 'r') as f:
        data = []
        datastream = f.read()
        datalines = [line.strip().split()
                     for line in datastream.split('\n') if line.strip() and line.strip().split()[0] != "#"]

        comment_index = []
        for line in datalines:
            try:
                comment_index.append(line.index('#'))
            except ValueError:
                comment_index.append(len(line))
                
        datalines = [line[0:comment_index[i]] for i, line in enumerate(datalines)]

        flags = list(zip(*datalines)[0])
        values = np.array([map(float, line) for line in zip(*datalines)[1:]]).T

        n_rows = len(values[0])
        n_data = len(values[:,0])
        if n_rows == 3:
            data = values
            data_covariances = np.array([[[0.] * 3] * 3] * n_data)
            data_covariances[:,0,0] = 1.
        elif n_rows == 6:
            data = values[:,0:3]
            data_covariances = np.array([[[0.] * 3] * 3] * n_data)
            data_covariances[:,0,0] = values[:,3]*values[:,3]
            data_covariances[:,1,1] = values[:,4]*values[:,4]
            data_covariances[:,2,2] = values[:,5]*values[:,5]
        elif n_rows == 9:
            data = values[:,0:3]
            data_covariances = np.array([[[0.] * 3] * 3] * n_data)
            data_covariances[:,0,0] = values[:,3]
            data_covariances[:,1,1] = values[:,4]
            data_covariances[:,2,2] = values[:,5]
            data_covariances[:,0,1] = values[:,6]
            data_covariances[:,1,0] = values[:,6]
            data_covariances[:,0,2] = values[:,7]
            data_covariances[:,2,0] = values[:,7]
            data_covariances[:,1,2] = values[:,8]
            data_covariances[:,2,1] = values[:,8]
        else:
            raise Exception('Your input file must have 4, 7, or 10 rows, where the first row is a string')


    return flags, data, data_covariances

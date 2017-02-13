#!/home/hamish/anaconda3/bin/python

import h5py 
import numpy as np
import sys
import os
print("HELLO")

direc = "InitialConditions"

if not os.path.exists(direc):
    os.makedirs(direc)

in_file = h5py.File("DATA/data.h5", 'r')

sg1 = "displacement"
sg2 = "east velocity"
sg3 = "north velocity"


data_eta = np.array(in_file[sg1][-1])
data_u = np.array(in_file[sg2][-1])
data_v = np.array(in_file[sg3][-1])

in_file.close()

init_file = h5py.File(direc+"/initial_conditions.h5", "w")

init_file.create_dataset(sg1, data=data_eta)
init_file.create_dataset(sg2, data=data_u)
init_file.create_dataset(sg3, data=data_v)

init_file.close()

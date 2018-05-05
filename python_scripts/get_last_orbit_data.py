import h5py
import numpy as np
import sys

load_file = str(sys.argv[1])

try:
    in_file = h5py.File(load_file, 'r')

except OSError:
    print("Cannot find data file " + load_file)
    sys.exit()

#try:
keys = list(in_file.keys())

data = []
for key in keys:
    data.append(in_file[key][-101:-2])

# Create new file
out_file = h5py.File(load_file.replace(".h5", "_orbit.h5"), 'w')

i = 0
for key in keys:
    out_file.create_dataset(key, data=data[i])
    i += 1

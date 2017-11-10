# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import numpy as np
import math
import h5py
from mpl_toolkits.basemap import Basemap
import sys

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Avante Garde'
plt.rcParams['mathtext.it'] = 'Avante Garde:italic'
plt.rcParams['mathtext.bf'] = 'Avante Garde:bold'
# plt.rc('font',sans_serif="Avante Garde")
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avante Garde']})
# plt.rc('font',serif='Palatino')
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)
plt.rc('text.latex',preamble='\\usepackage{siunitx},\\usepackage{sfmath}')

plt.rc('lines', linewidth=0.6)

plt.rc('figure', dpi=120)

# Creating a Triangulation without specifying the triangles results in the
# Delaunay triangulation of the points.
try:
    grid = np.loadtxt("input_files/grid_l6.txt",skiprows=1,usecols=(1,2))
except:
    d = input("Grid mismatch! Which grid do you want?: ")
    grid = np.loadtxt("input_files/grid_l" + str(int(d)) + ".txt",skiprows=1,usecols=(1,2))

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])

n = int(sys.argv[1])

in_file = h5py.File("DATA/data.h5", 'r')
# in_file = h5py.File("DATA_G6_T1000/data.h5", 'r')

data_u = np.array(in_file["east velocity"][:,0])
data_v = np.array(in_file["north velocity"][:,0])
data_eta = np.array(in_file["displacement"][:,0])

data_diss = np.array(in_file["avg dissipation output"][:])
plt.semilogy(data_diss)
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15,4))

t_max = 100
t_min = 0
time = np.linspace(t_min, t_max, len(data_eta))

# print(time)

ax1.plot(time, data_u)
ax2.plot(time, data_v)
ax3.plot(time, data_eta)
#
data_u = np.array(in_file["east velocity"][n])
data_v = np.array(in_file["north velocity"][n])
data_eta = np.array(in_file["displacement"][n])
data_diss = np.array(in_file["dissipated energy"][n])

print("TRIANGULATING POSITIONS")
# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, figsize=(15,4))
ax1.tricontourf(triang,data_u,11)
ax2.tricontourf(triang,data_v,11)
ax3.tricontourf(triang,data_eta,11)


fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)

m = Basemap(projection='npstere',boundinglat=80,lon_0=0)
x, y = m(np.degrees(x), np.degrees(y))

m.contourf(x, y, data_u, tri=True,ax=ax1)
m.contourf(x, y, data_v, tri=True,ax=ax2)
m.contourf(x, y, data_eta, tri=True,ax=ax3)


plt.show()

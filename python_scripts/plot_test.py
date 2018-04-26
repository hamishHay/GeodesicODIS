import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import numpy as np
import math
import h5py
from mpl_toolkits.basemap import Basemap
import sys
import math

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

n = int(sys.argv[1])

in_file = h5py.File("DATA/data.h5", 'r')
# in_file = h5py.File("DATA_G6_T1000/data.h5", 'r')

data_u = in_file["east velocity"][n]
data_v = in_file["north velocity"][n]

data_eta = in_file["displacement"][n]

N = int(1 + 0.5 * math.log((len(data_u)-2)/10, 2))

try:
    grid = np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))
except:
    d = input("Can't find grid file for level " + str(N) +"!!!")
    sys.exit()

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])

# data_diss = in_file["dissipated energy"][n]

print(np.max(data_u), np.max(data_v))#, np.max(data_eta))
print("TRIANGULATING POSITIONS")
# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, figsize=(15,4))
c1 = ax1.tricontourf(triang,data_u,21)
c2 = ax2.tricontourf(triang,data_v,21)
c3 = ax3.tricontourf(triang,data_eta,21)

plt.colorbar(c1, ax = ax1, orientation="horizontal")
plt.colorbar(c2, ax = ax2, orientation="horizontal")
plt.colorbar(c3, ax = ax3, orientation="horizontal")

fig.savefig("/home/hamish/Dropbox/plottings.pdf")
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)

m = Basemap(projection='npstere',boundinglat=80,lon_0=0)
x, y = m(np.degrees(x), np.degrees(y))

m.contourf(x, y, data_u, tri=True,ax=ax1)
m.contourf(x, y, data_v, tri=True,ax=ax2)
# m.contourf(x, y, data_eta, tri=True,ax=ax3)


plt.show()

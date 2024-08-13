import numpy as np
import h5py
import os
import sys
# from scipy.integrate import simps

import matplotlib as mpl
if "SSH_CONNECTION" in os.environ:
    mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1 import ImageGrid

# plt.rc('lines', linewidth=2.0)
#
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)



plt.rc('axes', titlesize=8)
plt.rc('axes', labelsize=7)
# plt.rc('axes', labelpad=10.0)
# plt.rc('axes', linewidth=1.20)



start = 675
end   = start+1
N=8

in_file = h5py.File("./DATA/data.h5", 'r')

d_u = in_file["east velocity"][start:end]
d_v = in_file["north velocity"][start:end]
d_h = in_file["displacement"][start:end]
f_lons = in_file["face longitude"][:]
f_lats = in_file["face latitude"][:]
grid = np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))
n_lons = grid[:,1]
n_lats = grid[:,0]

tri_face = tri.Triangulation(f_lons, f_lats)
tri_node = tri.Triangulation(n_lons, n_lats)

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(11,3.5))
axes = [ax1, ax2, ax3]


c1 = ax1.tricontourf(tri_face, d_u[0]/1e3, 11)
ax1.tricontour(tri_face, d_u[0]/1e3, 11, colors="k", linewidths=0.4)
cb1 = plt.colorbar(c1, orientation="horizontal", label="Eastward velocity [km/s]")
cb1.set_label(label="Eastward velocity [km/s]", size=8)

c2 = ax2.tricontourf(tri_face, d_v[0]/1e3, 11)
ax2.tricontour(tri_face, d_v[0]/1e3, 11, colors="k", linewidths=0.4)
cb2 = plt.colorbar(c2, orientation="horizontal")
cb2.set_label(label="Northward velocity [km/s]", size=8)

c3 = ax3.tricontourf(tri_node, d_h[0]/1e3, levels=np.linspace(-10, 10, 21), cmap=plt.cm.coolwarm)
ax3.tricontour(tri_node, d_h[0]/1e3, levels=np.linspace(-10, 10, 21), colors="k", linewidths=0.4)
cb3 = plt.colorbar(c3, orientation="horizontal")
cb3.set_label(label="Tidal height [km]", size=8)

for ax in axes:
    ax.set_aspect("equal")
    ax.set_xlabel("Longitude [°]", fontsize=8)
    ax.set_ylabel("Latitude [°]", fontsize=8)

fig.suptitle("$h_0 = 4$ km, $Ω = 2.91×10^{-4}$ rad s$^{-1}$, $M_i = 6.42×10^{23}$ kg", fontsize=10,y=0.8)

fig.savefig("giant_impact_test.pdf",bbox_inches="tight")

plt.show()

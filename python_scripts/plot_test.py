import matplotlib as mpl
import os
if "SSH_CONNECTION" in os.environ:
    mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import numpy as np
import math
import h5py
# from mpl_toolkits.basemap import Basemap
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

# data_u = np.mean(in_file["east velocity"][-10*1000:-1], axis=0)
# data_v = np.mean(in_file["north velocity"][-10*1000:-1], axis=0)

# data_diss = in_file["east velocity"][:]**2 + in_file["north velocity"][:]**2
# data_diss = 1.2e3 * 1000 * np.mean(data_diss, axis=0)

# print(np.max(in_file["displacement"][-3*1000:]))
data_eta = in_file["displacement"][n]
# data_eta = data_diss
#data_eta = data_diss

# data_eta = np.mean(in_file["dissipated energy"][:], axis=0)

N = int(1 + 0.5 * math.log((len(data_u)-2)/10, 2))

try:
    grid = np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))
except:
    d = input("Can't find grid file for level " + str(N) +"!!!")
    sys.exit()

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])

# data_diss = in_file["dissipated energy"][n]

print(np.max(abs(data_u)), np.max(abs(data_v)), np.max(data_eta))
print("TRIANGULATING POSITIONS")
# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, figsize=(15,4))
c1 = ax1.tripcolor(triang,data_u)
# c2 = ax2.tripcolor(triang,data_v)
ax2.quiver(triang.x, triang.y, data_u, data_v)
c3 = ax3.tripcolor(triang,data_eta)
ax3.triplot(triang)
# c1 = ax1.tricontourf(triang,data_u,21)
# c2 = ax2.tricontourf(triang,data_v,21)
# c3 = ax3.tricontourf(triang,data_eta,21)

plt.colorbar(c1, ax = ax1, orientation="horizontal")
# plt.colorbar(c2, ax = ax2, orientation="horizontal")
plt.colorbar(c3, ax = ax3, orientation="horizontal")

fig.savefig("/home/hamish/Dropbox/Tests/plottings.pdf")

# import cartopy
# import cartopy.crs as ccrs
#
# from scipy.interpolate import SmoothSphereBivariateSpline
# def smooth_data(x, y, xx, yy, data_u, data_v,name=None,s=0.05):
#     print("Interpolating u velocity")
#
#     max_u, max_v = [np.amax(data_u), np.amax(data_v)]
#     interp_u = SmoothSphereBivariateSpline(np.radians(90.-y), np.radians(x), data_u/max_u, s=s)
#
#     print("Interpolating v velocity")
#     interp_v = SmoothSphereBivariateSpline(np.radians(90.-y), np.radians(x), data_v/max_v, s=0.05)
#
#     print("Generating interpolated data")
#     data_interp_u = interp_u(yy, xx)*max_u
#     data_interp_v = interp_v(yy, xx)*max_v
#     if name != None:
#         np.savetxt(name + "/data_u_interp.txt", data_interp_u)
#         np.savetxt(name + "/data_v_interp.txt", data_interp_v)
#     return data_interp_u, data_interp_v
#
#
# nx, ny = 301, 151
# xx, yy = [np.linspace(0, 360, nx), np.linspace(0, 180, ny)]
# xx, yy = np.radians(xx), np.radians(yy)
#
#
# # proj = ccrs.AzimuthalEquidistant(central_latitude=-90)
# # proj = ccrs.SouthPolarStereo()
# proj=ccrs.PlateCarree()
# ax = plt.axes(projection=proj)#projection=ccrs.PlateCarree())
# fig = plt.gcf()
# fig.set_size_inches(6, 3.5)
#
# data_u = in_file["east velocity"][-100*10:]
# data_v = in_file["north velocity"][-100*10:]
#
# data_eta = 1e3 * 4e-3 * np.sqrt(data_u**2. + data_v**2.)*(data_u**2. + data_v**2.)
# data_eta = np.mean(data_eta, axis=0)
#
# data_u = in_file["east velocity"][n]
# data_v = in_file["north velocity"][n]
#
# data_crs = ccrs.PlateCarree()
# du, dv = smooth_data(np.degrees(x), np.degrees(y), xx, yy, data_u, data_v)
# deta, dv = smooth_data(np.degrees(x), np.degrees(y), xx, yy, data_eta, data_v)
# # crs = ccrs.Orthographic()
# xx, yy = [np.linspace(-180, 180., nx), np.linspace(-90., 90., ny)]
# # ax.quiver(xx, yy, du, dv,  transform=data_crs, regrid_shape=30)
#
# deta[deta<0.0] = 0.0
# levels = np.linspace(np.amin(deta*1e3), np.amax(deta*1e3), 16)
# #levels = np.linspace(0.,4.2,11)
# #ticks = np.linspace(1, 4., 7)
# c1 = ax.contourf(xx, yy, deta*1e3, levels=levels, transform=data_crs, cmap='magma')
# cb = plt.colorbar(c1, ax = ax, orientation="vertical", shrink=0.68)
# cb.set_label("Heat flux [$\\num{e-8} \si{\milli\watt\per\square\metre}$]")
#
# # ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
# ax.set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
# ax.set_yticks([-90, -45, 0, 45, 90])
# ax.set_xlabel("Longitude [\si{\degree}]")
# ax.set_ylabel("Latitude [\si{\degree}]")
#
# ax.set_title("Eccentricity and obliquity forcing")
#
# fig.savefig("heat_flux_map_obliq.pdf", bbox_inches='tight')



plt.show()

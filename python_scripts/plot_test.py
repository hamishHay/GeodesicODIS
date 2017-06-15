import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import h5py
from mpl_toolkits.basemap import Basemap
import sys

def potential_ecc_dlat(lat,lon,radius,omega,e,time):
    # Convert lat to colat
    # lat = np.pi/2.0 - lat

    # If numpy array, convert lat and lon using meshgrid
    # lat, lon = np.meshgrid(lat,lon)

    factor = omega**2 * radius**2*e


    M = omega*time
    pot = -factor *0.75*np.sin(2*lat) * (3*np.cos(M) * (1.0 + np.cos(2*lon)) + 4*np.sin(M)*np.sin(2*lon))


    return pot

# Creating a Triangulation without specifying the triangles results in the
# Delaunay triangulation of the points.
grid = np.loadtxt("input_files/grid_l6.txt",skiprows=1,usecols=(1,2))

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])

n = int(sys.argv[1])

in_file = h5py.File("DATA/data.h5", 'r')

data_u = np.array(in_file["east velocity"][n])
data_v = np.array(in_file["north velocity"][n])

print(np.amax(data_u), np.amax(data_v))

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(10,4))
#
# lats = y
# lons = x
# data = potential_ecc_dlat(lats,lons,2575.5e3,4.559e-6,0.0288, 0.0)

# map = Basemap(projection='hammer',lon_0=0,lat_0 = 0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
# # map = Basemap(projection='hammer',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
# x, y = map(np.rad2deg(x), np.rad2deg(y))
# cs = map.contourf(x,y,data_u,31,cmap=plt.cm.Spectral, tri=True)#,vmax=2400)
# cbar = map.colorbar(cs,location='bottom',pad="5%")
#
# map = Basemap(projection='hammer',lon_0=0,lat_0 = 0, ax=ax2)#, llcrnrlon=-90,urcrnrlon=90)
# cs = map.contourf(x,y,data_v,31,cmap=plt.cm.Spectral, tri=True)#,vmax=2400)
# cbar = map.colorbar(cs,location='bottom',pad="5%")


tcnt = ax1.tricontourf(triang, data_v, 21)
cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')

tcnt = ax2.tricontourf(triang, data_u, 21)
cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')


ax1.set_aspect('equal')
ax2.set_aspect('equal')

plt.show()

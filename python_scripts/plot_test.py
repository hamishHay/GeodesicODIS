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
try:
    grid = np.loadtxt("input_files/grid_l6.txt",skiprows=1,usecols=(1,2))
except:
    d = input("Grid mismatch! Which grid do you want?: ")
    grid = np.loadtxt("input_files/grid_l" + str(int(d)) + ".txt",skiprows=1,usecols=(1,2))

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])

n = int(sys.argv[1])

in_file = h5py.File("DATA/data.h5", 'r')

data_u = np.array(in_file["east velocity"][n])
data_v = np.array(in_file["north velocity"][n])
data_eta = np.array(in_file["displacement"][n])
# data_diss = np.array(in_file["dissipated energy"][n])
data_diss = np.mean(np.array(in_file["dissipated energy"]),axis=0)
# data_u = -np.cos(y)/252.1e3 * 0.11

data_t_diss = np.array(in_file["avg dissipation output"])[:-1] * 4* np.pi * (2575.5e3)**2



print("u_max: "+ str(np.amax(abs(data_u))), "v_max: " + str(np.amax(abs(data_v))), "eta_max: " + str(np.amax(data_eta)), "diss_avg: " + str(np.mean(data_t_diss[-101:])))

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, ax1 = plt.subplots()
time = np.linspace(0, 40, len(data_t_diss))
ax1.semilogy(time,data_t_diss,lw=0.5)
# ax1.semilogy(time,np.mean(data_diss))

fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=1,ncols=4,figsize=(12,3),dpi=120)

tcnt = ax1.tricontourf(triang, data_v, 21, cmap = plt.cm.winter)
cnt = ax1.tricontour(triang, data_v, 11, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')

tcnt = ax2.tricontourf(triang, data_u, 21, cmap = plt.cm.winter)
cnt = ax2.tricontour(triang, data_u, 11, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')

tcnt = ax3.tricontourf(triang, data_eta, 21, cmap = plt.cm.winter)
cnt = ax3.tricontour(triang, data_eta, 11, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')

tcnt = ax4.tricontourf(triang, data_diss, 21, cmap = plt.cm.plasma)
cnt = ax4.tricontour(triang, data_diss, 11, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax4, orientation='horizontal')

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

ax1.set_xlim([0,2*np.pi])
ax2.set_xlim([0,2*np.pi])
ax3.set_xlim([0,2*np.pi])
ax4.set_xlim([0,2*np.pi])

ax1.set_ylim([-np.pi*0.5,np.pi*0.5])
ax2.set_ylim([-np.pi*0.5,np.pi*0.5])
ax3.set_ylim([-np.pi*0.5,np.pi*0.5])
ax4.set_ylim([-np.pi*0.5,np.pi*0.5])


ax1.set_title("North Vel")
ax2.set_title("East Vel")
ax3.set_title("Displacement")
ax4.set_title("Dissipated Energy")

plt.tight_layout()
# plt.show()


fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=1,ncols=5,figsize=(15,3),dpi=120)

lat0 = 20
lon0 = 180
proj = 'ortho'

x = np.degrees(x)
y = np.degrees(y)

parallels = np.arange(0.,81,20.)
meridians = np.arange(10.,351.,30.)

# map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
# # map = Basemap(projection=proj,llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
x, y = map(x, y)
cs = map.contourf(x,y,data_v,31,cmap=plt.cm.winter, tri=True)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels)
map.drawmeridians(meridians)

map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax2)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_u,31,cmap=plt.cm.winter, tri=True)
# map.quiver(x,y,data_u,data_v,scale=3)
# cnt = map.contour(x,y,data_u,levels=[0.0], tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels)
map.drawmeridians(meridians)

map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax3)#, llcrnrlon=-90,urcrnrlon=90)
map.quiver(x,y,data_u,data_v,scale=0.1)
map.drawparallels(parallels)
map.drawmeridians(meridians)
map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax4)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_eta,31,cmap=plt.cm.winter, tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels)
map.drawmeridians(meridians)

map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax5)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_diss,31,cmap=plt.cm.plasma, tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels)
map.drawmeridians(meridians)

plt.tight_layout()

fig, ax = plt.subplots()

x = np.radians(grid[:,1])
y = np.radians(grid[:,0])
# data_eta[data_eta == 0] = -5

ax.scatter(x,y,c=data_eta)
ax.set_aspect("equal")

plt.show()

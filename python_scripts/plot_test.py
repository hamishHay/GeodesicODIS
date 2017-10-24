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
# in_file = h5py.File("DATA/data.h5", 'r')

data_u = np.array(in_file["east velocity"][n])
data_v = np.array(in_file["north velocity"][n])
data_eta = np.array(in_file["displacement"][n])
data_diss = np.array(in_file["dissipated energy"][n])


x = x[2:]
y = y[2:]
data_u = data_u[2:]
data_v = data_v[2:]
data_eta = data_eta[2:]
data_diss = data_diss[2:]



print("DATA LOADED")

r = 2.634100E+06
# r = 252.1e3
#
# u_test = np.cos(y)**2.0 / r * (np.sin(2*x) + np.sin(4*x))
#
# v_test = 3 * np.cos(2*x)*np.cos(x)**2.0 * np.cos(y)**2.0 * np.sin(y) / r
#
# data_eta = v_test #data_v
# data_diss = u_test #data_u
#
# # data_u = abs(data_u - u_test)
# # data_v = abs(data_v - v_test)
#
# data_u[np.isnan(data_u)] = 0.0
# data_v[np.isnan(data_v)] = 0.0


# data_diss = -np.cos(2*x)**2.0 * np.sin(2*y)/(r*np.cos(y)) - np.sin(x)/r
#
# data_diss = data_diss - data_eta
#
# data_eta *= 1000.0

# print(np.amax(data_diss[-100:]))

# data_diss = np.mean(np.array(in_file["dissipated energy"][-100*20:,2:]),axis=0)
# data_diss = np.mean(np.array(in_file["dissipated energy"][-100:,:]),axis=0)
# data_u = -np.cos(y)/252.1e3 * 0.11
# data_diss = np.mean(np.array(in_file["dissipated energy"][-100:,:])
data_t_diss = np.array(in_file["avg dissipation output"])[:-1] * 4* np.pi * r**2



# print("u_max: "+ str(np.amax(abs(data_u))), "v_max: " + str(np.amax(abs(data_v))), "eta_max: " + str(np.amax(data_eta)), "diss_avg: " + str(np.mean(data_t_diss[-101:])))

print("TRIANGULATING POSITIONS")
# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

t_max = len(data_t_diss)
t_min = 0

fig, ax1 = plt.subplots()
time = np.linspace(0, t_max, len(data_t_diss))
ax1.semilogy(time,data_t_diss,lw=1.0)
ax1.set_xlim([t_min,t_max])
ax1.grid()
ax1.set_title("Dissipation")


data_u2 = np.amax(np.array(in_file["east velocity"]),axis=1)
fig, ax1 = plt.subplots()
time = np.linspace(0, t_max, len(data_u2))
ax1.semilogy(time,data_u2,lw=1.0)
ax1.set_xlim([t_min,t_max])
ax1.grid()
ax1.set_title('Max vel')
# ax1.semilogy(time,np.mean(data_diss))

# fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=1,ncols=4,figsize=(12,3),dpi=120)
#
# fig, ax4 = plt.subplots(dpi=120)
#
# print("PLOTTING VELOCITY")
# # tcnt = ax1.tricontourf(triang, data_v, 21, cmap = plt.cm.winter)
# # cnt = ax1.tricontour(triang, data_v, 11, colors='k', linewidths=0.5)
# # cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
#
# # tcnt = ax2.tricontourf(triang, data_u, 21, cmap = plt.cm.winter)
# # cnt = ax2.tricontour(triang, data_u, 11, colors='k', linewidths=0.5)
# # cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
#
# print("PLOTTING PRESSURE")
# # tcnt = ax3.tricontourf(triang, data_eta, 21, cmap = plt.cm.winter)
# # cnt = ax3.tricontour(triang, data_eta, 11, colors='k', linewidths=0.5)
# # cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
#
# tcnt = ax4.tricontourf(triang, data_diss, 21, cmap = plt.cm.plasma)
# cnt = ax4.tricontour(triang, data_diss, 11, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax4, orientation='horizontal')
#
# # ax1.set_aspect('equal')
# # ax2.set_aspect('equal')
# # ax3.set_aspect('equal')
# ax4.set_aspect('equal')
#
# # ax1.set_xlim([0,2*np.pi])
# # ax2.set_xlim([0,2*np.pi])
# # ax3.set_xlim([0,2*np.pi])
# ax4.set_xlim([0,2*np.pi])
#
# # ax1.set_ylim([-np.pi*0.5,np.pi*0.5])
# # ax2.set_ylim([-np.pi*0.5,np.pi*0.5])
# # ax3.set_ylim([-np.pi*0.5,np.pi*0.5])
# ax4.set_ylim([-np.pi*0.5,np.pi*0.5])
#
#
# # ax1.set_title("North Vel")
# # ax2.set_title("East Vel")
# # ax3.set_title("Displacement")
# ax4.set_title("Dissipated Energy [W]")
#
# # ax1.set_title("North Vel Diff")
# # ax2.set_title("East Vel Vel Diff")
# # ax3.set_title("North Vel")
# # ax4.set_title("East Vel")
#
# plt.tight_layout()
# plt.show()
data_eta2 = np.array(in_file["displacement"][:])
max_p = np.amax(data_eta2,axis=1)

fig, ax = plt.subplots()

ax.plot(max_p)
ax.grid(which='both')


fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=1,ncols=5,figsize=(15,3),dpi=120)

lat0 = 0
lon0 = 0

proj = 'ortho'

x = np.degrees(x)
y = np.degrees(y)

parallels = np.arange(-80,81,20.)
meridians = np.arange(0,351.,30.)

lw=0.2

# map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
map = Basemap(projection=proj,boundinglat=60,lon_0=lon0,lat_0=lat0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
# # map = Basemap(projection=proj,llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
x, y = map(x, y)
cs = map.contourf(x,y,data_v,31,cmap=plt.cm.winter, tri=True)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels, linewidth=lw)
map.drawmeridians(meridians, linewidth=lw)

map = Basemap(projection=proj,boundinglat=60,lon_0=lon0,lat_0=lat0, ax=ax2)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_u,31,cmap=plt.cm.winter, tri=True)
# map.quiver(x,y,data_u,data_v,scale=3)
# cnt = map.contour(x,y,data_u,levels=[0.0], tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels, linewidth=lw)
map.drawmeridians(meridians, linewidth=lw)

map = Basemap(projection=proj,boundinglat=60,lon_0=lon0,lat_0=lat0, ax=ax3)#, llcrnrlon=-90,urcrnrlon=90)
map.quiver(x,y,data_u,data_v,scale=3e-4)
map.drawparallels(parallels, linewidth=lw)
map.drawmeridians(meridians, linewidth=lw)

map = Basemap(projection=proj,boundinglat=60,lon_0=lon0,lat_0=lat0, ax=ax4)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_eta,31,cmap=plt.cm.winter, tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels, linewidth=lw)
map.drawmeridians(meridians, linewidth=lw)

map = Basemap(projection=proj,lon_0=lon0,lat_0=lat0, ax=ax5)#, llcrnrlon=-90,urcrnrlon=90)
cs = map.contourf(x,y,data_diss,31,cmap=plt.cm.coolwarm, tri=True)#,vmax=2400)
cbar = map.colorbar(cs,location='bottom',pad="5%")
map.drawparallels(parallels, linewidth=lw)
map.drawmeridians(meridians,linewidth=lw)

plt.tight_layout()

# fig, ax1 = plt.subplots(dpi=120)
#
# data_p = np.array(in_file["displacement"][-101:-1]) + 1000.0 * 0.11 * 10e3
# # data_p = np.mean(data_p, axis=0)
#
# x = grid[:,1]
# y = grid[:,0]
#
# data_p = data_p[:, y <= -60.0]
#
# data_p = np.mean(data_p/1e6, axis=1)
#
# time = np.linspace(0,1,100)
#
# ax1.plot(time, data_p)
# ax1.set_title("Enceladus SPT Averaged Pressure")
# ax1.set_ylabel("Pressure [\si{\mega\pascal}]")
# ax1.set_xlabel("Time $t/T$")
#
# ax1.grid(which='major', alpha=0.8)
# ax1.grid(which='minor', alpha=0.2)


# fig, ax1 = plt.subplots(dpi=120)
#
#
#
# proj = 'spstere'
# # grid2 = np.loadtxt("input_files/grid_l7.txt",skiprows=1,usecols=(1,2))
#
# x = grid[:,1]
# y = grid[:,0]
#
# x = x[2:]
# y = y[2:]
#
# map = Basemap(projection=proj,boundinglat=-60,lon_0=lon0,lat_0=lat0, ax=ax1)#, llcrnrlon=-90,urcrnrlon=90)
# x, y = map(x, y)
#
# cs = map.contourf(x,y,data_diss,12,cmap=plt.cm.plasma, tri=True)
# cbar = map.colorbar(cs,location='bottom',pad="5%",label='Dissipation [W]')
# #
# # parallels = np.arange(-80,81,10.)
# #
# # dmin = -np.amax(abs(data_eta))
# # dmax = -dmin
# # norm = mpl.colors.Normalize(vmin=dmin,vmax=dmax)
# #
# # avg_p = np.mean(np.array(in_file["displacement"]), axis=0)
# #
# # # cs = map.contour(x,y,data_eta[:len(x)],31,colors='k', tri=True)
# # map.quiver(x,y,data_u[:len(x)],data_v[:len(x)],scale=4)
# map.drawparallels(parallels, linewidth=lw)
# map.drawmeridians(meridians,labels=[True,True,True,True],linewidth=lw)


# fig.savefig("obliq_cd_polar.pdf",bbox_inches='tight')
# fig, ax = plt.subplots()
#
# x = np.radians(grid[:,1])
# y = np.radians(grid[:,0])
# # data_eta[data_eta == 0] = -5
#
# ax.scatter(x,y,c=data_eta)
# ax.set_aspect("equal")

plt.show()

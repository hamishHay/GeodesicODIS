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
from pyshtools.expand import SHExpandLSQ
import pyshtools
import cartopy.crs as ccrs
from cartopy.vector_transform import vector_scalar_to_grid


plt.style.use("dark_background")
plt.rcParams['axes.facecolor'] = 'k'
plt.rcParams['contour.negative_linestyle']= 'solid'

# plt.rcParams['mathtext.fontset'] = 'custom'
# plt.rcParams['mathtext.rm'] = 'Avante Garde'
# plt.rcParams['mathtext.it'] = 'Avante Garde:italic'
# plt.rcParams['mathtext.bf'] = 'Avante Garde:bold'
# # plt.rc('font',sans_serif="Avante Garde")
# plt.rc('font',**{'family':'sans-serif','sans-serif':['Avante Garde']})
# # plt.rc('font',serif='Palatino')
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rc('text', usetex=True)
# plt.rc('text.latex',preamble='\\usepackage{siunitx},\\usepackage{sfmath}')

plt.rc('lines', linewidth=0.6)

plt.rc('figure', dpi=100)

# Creating a Triangulation without specifying the triangles results in the
# Delaunay triangulation of the points.

nn = int(sys.argv[1])
in_file = h5py.File("DATA/data.h5", 'r')

for i in range(nn, nn+1):
    n = i
    # in_file = h5py.File("DATA_G6_T1000/data.h5", 'r')

    data_u = in_file["east velocity"][n]
    data_v = in_file["north velocity"][n]
    vel_lon = in_file["face longitude"][:] - 180.0
    vel_lat = in_file["face latitude"][:]

    data_eta = in_file["displacement"][n]

    input_coord_system = ccrs.PlateCarree()

    this_plot_projection =  ccrs.PlateCarree()

    (v_lon, v_lat, v_u, v_v, v_speed) = vector_scalar_to_grid(input_coord_system,     \
                            this_plot_projection, \
                            (20,10), \
                            vel_lon, vel_lat, \
                            data_u, \
                            data_v, \
                            data_v )
    #
    # print(vel_lon)

    # data_u = np.mean(in_file["east velocity"][-10*1000:-1], axis=0)
    # data_v = np.mean(in_file["north velocity"][-10*1000:-1], axis=0)

    # data_diss = in_file["east velocity"][:]**2 + in_file["north velocity"][:]**2
    # data_diss = 1.2e3 * 1000 * np.mean(data_diss, axis=0)
    
    # print(np.max(in_file["displacement"][-3*1000:]))
    N = int(1 + 0.5 * math.log((len(data_eta)-2)/10, 2))
    data_eta = data_eta #- np.mean(in_file["displacement"][-2001:-2], axis=0)
    lmin = np.amin(in_file["displacement"][n:n+100,:])
    lmax = np.amax(in_file["displacement"][n:n+100,:])
    norm = max(abs(lmin), abs(lmax))
    
    levels = np.linspace(lmin, lmax, 13)

    P = 2*np.pi/(2.05e-5)

    # data_u = np.mean(in_file["east velocity"][-101:-2], axis=0)
    # data_diss = np.mean(in_file["dissipated energy"][-101:-2], axis=0)


    try:
        grid = np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))
    except:
        d = input("Can't find grid file for level " + str(N) +"!!!")
        sys.exit()

    x = grid[:,1] - 180.0
    y = grid[:,0]

    colat = np.radians(y)
    lons = np.radians(y)
   
    mag = np.amax(abs(np.sqrt(in_file["east velocity"][n:n+100]**2.0 + in_file["north velocity"][n:n+100]**2.0)))
    # reg = mag > 1e-7
    # # print(reg)
    # x = x[reg]
    # y = y[reg]
    # # for i in range(len(reg)):
    # #     print(x[i], y[i])
    #
    # # print(x, y)
    # data_eta = data_eta[reg]
    # data_u = data_u[reg]
    # data_v = data_v[reg]
    # region = region[region]

    # data_diss = in_file["dissipated energy"][n]

    # print(np.max(abs(data_u)), np.max(abs(data_v)), np.amax(data_eta))
    # print(np.max(abs(data_u)), np.max(abs(data_v)), np.amin(data_eta))
    print("TRIANGULATING POSITIONS")
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y)
    # for i,j in zip(vel_lon, vel_lat):
    #     print(i,j)
    tri_vel = tri.Triangulation(vel_lon, vel_lat)

    

    fig, ax1 = plt.subplots(ncols=1, figsize=(6,4))#,subplot_kw={'projection': ccrs.PlateCarree()})
    
    
    # c2 = ax2.tricontourf(tri_vel,data_diss,11,transform=ccrs.PlateCarree())
    print(lmin, lmax)
    norm = mpl.colors.Normalize(vmin=lmin, vmax=lmax)
    
    c1 = [ax1.tricontourf(triang,data_eta,levels)]#, norm=norm, transform=ccrs.PlateCarree())]
    # Q = ax1.quiver(tri_vel.x, tri_vel.y, data_u/mag, data_v/mag, scale=40.0,transform=ccrs.PlateCarree(),regrid_shape=25, pivot='mid')
    # c1 = ax1.scatter(x, y, c=data_eta)
    # c2 = [ax1.tricontour(triang, data_eta, levels=levels, colors=[(0.2, 0.2, 0.2)])]
    plt.colorbar(c1[0], ax = ax1, orientation="horizontal", label='Tidal Amplitude [m]')
    # c2 = 
    # gridlines = ax1.gridlines(draw_labels=True)
    # Q = ax1.quiver(v_lon, v_lat, v_u, v_v,zorder=10, scale=1.0, color=(0.1,0.1,0.1))

    # gridlines.xlabels_top = False
    # gridlines.ylabels_right = False
    # gridlines.xlines = False
    # gridlines.ylines = False
    # gridlines.xlocator = mpl.ticker.FixedLocator(np.arange(-180,181,45))

    # ax1.set_facecolor("k")
    # ax1.background_patch.set_facecolor('k')
    # ax1.outline_patch.set_edgecolor('w')
from potential import potential_ecc 
# LON, COLAT = np.meshgrid(lons, colat)    
times = in_file['kinetic avg output']

def update(i):
    print(i)
    data_new = in_file["displacement"][n+i]
    data_u = in_file["east velocity"][n+i]
    data_v = in_file["north velocity"][n+i]

    (v_lon, v_lat, v_u, v_v, v_speed) = vector_scalar_to_grid(input_coord_system,     \
                            this_plot_projection, \
                            (20,10), \
                            vel_lon, vel_lat, \
                            data_u, \
                            data_v, \
                            data_v )
   
    time = times[n+i] 
    # print(time)
    U = potential_ecc(0.0094, 1565.e3, 2.05e-5, time, np.radians(90-triang.y), np.radians(triang.x))


    # ax1.set_facecolor("k")
    for tp in c1[0].collections:
        tp.remove()
    c1[0] = ax1.tricontourf(triang,data_new,levels=levels)#transform=ccrs.PlateCarree())
    ax1.set_title("{:1.1f} Days".format(i/100.0 * 3.55))

    # for tp in c2[0].collections:
    #     tp.remove()

    # c2[0] = ax1.tricontour(triang, U/1.3, levels=levels, colors=[(0.05, 0.05, 0.05)])
    Q = ax1.quiver(v_lon, v_lat, v_u, v_v)

    Q.set_UVC(v_u, v_v)

    ax1.set_xlim([-180,180])
    ax1.set_ylim([-90,90])
    # Q = ax1.quiver(tri_vel.x, tri_vel.y, data_u/mag, data_v/mag, scale=40.0,transform=ccrs.PlateCarree(),regrid_shape=17, pivot='mid', color=(0.1,0.1,0.1))
    return  c1[0].collections, Q#, c2[0].collections


from matplotlib import animation
ani = animation.FuncAnimation(fig, update, 100, interval=60)

ani.save("/home/hamish/Dropbox/Tests/plottings.gif", writer='imagemagick', dpi=120)

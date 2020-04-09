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
    #
    # print(vel_lon)

    # data_u = np.mean(in_file["east velocity"][-10*1000:-1], axis=0)
    # data_v = np.mean(in_file["north velocity"][-10*1000:-1], axis=0)

    data_diss = np.mean(in_file["east velocity"][-101:-2]**2 + in_file["north velocity"][-101:-2]**2, axis=0)
    # data_diss = 1.2e3 * 1000 * np.mean(data_diss, axis=0)
    data_eta = in_file["displacement"][n]
    # print(np.max(in_file["displacement"][-3*1000:]))
    N = int(1 + 0.5 * math.log((len(data_eta)-2)/10, 2))
    data_eta = data_eta #- np.mean(in_file["displacement"][-2001:-2], axis=0)
    levels=np.linspace(10, 380,21)#np.amax(data_eta), 21)
    # d = np.ma.array(d, mask=d < .2)
    # data_eta = np.ma.array(data_eta, mask=data_eta < 10)

    P = 2*np.pi/(2.05e-5)
    # print(in_file["kinetic avg output"][n]/P)
    # data_eta = data_diss
    #data_eta = data_diss

    # data_u = np.mean(in_file["east velocity"][-101:-2], axis=0)
    # data_diss = np.mean(in_file["dissipated energy"][-101:-2], axis=0)


    try:
        grid = np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))
    except:
        d = input("Can't find grid file for level " + str(N) +"!!!")
        sys.exit()

    x = grid[:,1] - 180.0
    y = grid[:,0]


    # region = np.ones(len(x), dtype=np.int)
    # # region[x<270] = 0
    # region[(x>120) & (x<270)] = 0
    # # # region[y<80] = 0
    # # # region[y>-80] = 0
    # reg = (x>180) & (x<270) & (y<80) & (y>-80)
    #
    mag = abs(data_u**2.0 + data_v**2.0)
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
    print(len(data_diss), len(data_u))
    # x, y, u, v, c = vector_scalar_to_grid(ccrs.Mollweide(), ccrs.Mollweide(), 40, x, y, data_eta, data_eta,)

    # mask = np.zeros(len(x), np.int)
    # mask[137] = 1
    # mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
    # triang.set_mask(mask)

    # fig, ax = plt.subplots()
    # levels = np.arange(-0.20, 0.32+0.01, 0.04)
    # levels = np.arange(-0.06, 0.1+0.01, 0.02)
    # levels = np.arange(-0.2, 0.32+0.01, 0.04)
    # c = ax.tricontourf(triang,data_eta*1e6)#, colors='k', ls='-')
    # plt.colorbar(c, ax = ax, orientation="horizontal", label='Heat flux [\\num{e-6} \\si{\\watt\\per\\metre\\squared}]')
    # ax.set_xlim([0, 360])
    # ax.set_ylim([-90, 90])
    # ax.set_aspect('equal')
    #
    # fig.savefig("/home/hamish/figure_odis.png", bbox_inches='tight')

    # plt.colorbar(c, ax = ax, orientation="horizontal")

    # cilm, chi2 = SHExpandLSQ (data_eta, y, x, 12, norm=1)
    #
    # print(cilm[0])
    #
    # power_per_l = pyshtools.spectralanalysis.spectrum(cilm)
    # degrees = np.arange(cilm.shape[1])
    #
    # fig, ax = plt.subplots(1, 1)
    # ax.plot(degrees, power_per_l, 'o-')
    # ax.set(yscale='log', xscale='log', xlabel='Spherical harmonic degree', ylabel='Power')
    # # ax.grid()
    # fig.savefig("/home/hamish/Dropbox/Tests/eta_spectrum.pdf", bbox_inches='tight')

    # ax.plot(np.sum(cilm[0]**2.0, axis=-1), 'o')

    fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, figsize=(15,4),subplot_kw={'projection': ccrs.PlateCarree()})
    # c1 = ax1.tripcolor(triang,data_u)
    # c2 = ax2.tripcolor(triang,data_v)
    # c3 = ax3.tripcolor(triang,data_eta)
    # ax3.triplot(triang)
    # ax3.patch.set_color('.25')
    # c1 = ax1.tricontourf(triang,data_u,11)
    # c2 = ax2.tricontourf(triang,data_v,11)
    # print(len(triang), len(data_diss))
    c1 = ax1.tricontourf(tri_vel,data_u,11,transform=ccrs.PlateCarree())
    c2 = ax2.tricontourf(tri_vel,data_v,11,transform=ccrs.PlateCarree())
    c3 = ax3.tricontourf(tri_vel,data_diss,11,transform=ccrs.PlateCarree())
    # ax3.quiver(tri_vel.x, tri_vel.y, data_u/mag, data_v/mag, scale=40.0,transform=ccrs.PlateCarree(),regrid_shape=25, pivot='mid')
    # c3 = ax3.scatter(x, y, c=data_eta)

    plt.colorbar(c1, ax = ax1, orientation="horizontal")
    plt.colorbar(c2, ax = ax2, orientation="horizontal")
    plt.colorbar(c3, ax = ax3, orientation="horizontal")

    # ax2.set_xlim([160, 200])
    # ax2.set_ylim([-30, 30])
    # ax3.set_xlim([90, 270])
    # ax3.set_ylim([-90, 90])

    fig.savefig("/home/hamish/Dropbox/Tests/plottings.pdf")
    # fig.savefig("/home/hamish/Dropbox/Tests/plottings_{:03d}.png".format(n), bbox_inches='tight', dpi=300)

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

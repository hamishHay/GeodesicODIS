import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import h5py
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from scipy.interpolate import griddata
# import cmocean

plt.rc('font', family='serif')
plt.rc('font',serif='Palatino')
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)
plt.rc('text.latex',preamble='\\usepackage{siunitx}')

# Set matplotlib rc options
plt.rcParams['axes.linewidth'] = 0.8



n = int(sys.argv[1])
N = int(sys.argv[2])
nx = 201
ny = int(nx/2)

grid = np.loadtxt("input_files/grid_l" +str(N)+ ".txt",skiprows=1,usecols=(1,2))


x = np.radians(grid[:,1])
y = np.radians(grid[:,0])
points = np.fliplr(np.radians(grid))

lons = np.linspace(0, np.pi*2, nx)
lats = np.linspace(np.pi*0.5, -np.pi * 0.5, ny)

# ax.xaxis.set_ticks(np.arange(start, end, stepsize))
xticks = np.arange(0,360,45)
yticks = np.arange(-90,91,45)

points2 = []
for i in range(len(lons)):
    for j in range(len(lats)):
        points2.append([lons[i], lats[j]])
points2 = np.array(points2)
# in_file = h5py.File("DATA_H5000_ECC/data.h5", 'r')
# in_file = h5py.File("DATA_H5000_OBLIQ/data.h5", 'r')
fig, axes1 = plt.subplots(nrows=2,ncols=4,sharex='col',sharey='row', figsize=(15,5),dpi=120)

# axes = [axes1[0][0], axes1[0][1], axes1[1][0], axes1[1][1]]
axes = [axes1[0][0], axes1[0][1], axes1[0][2], axes1[0][3], axes1[1][0], axes1[1][1], axes1[1][2], axes1[1][3]]

l_axes = [axes1[0][0], axes1[1][0]]
b_axes = [axes1[1][0], axes1[1][1], axes1[1][2], axes1[1][3]]

time = np.array([0, 12, 25, 37, 50, 62, 75, 87], dtype=np.int) + n
time = np.array([0, 6, 12, 18, 25, 31, 37, 43], dtype=np.int) + n
times = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]

in_file = h5py.File("DATA/data.h5", 'r')

dmin = 0.0
dmax = 2.1e-5

norm = mpl.colors.Normalize(vmin=dmin,vmax=dmax)
levels = np.linspace(dmin,dmax,201)
for i in range(len(time)):

    data_u = np.array(in_file["east velocity"][time[i]])
    data_v = np.array(in_file["north velocity"][time[i]])

    interp_u = griddata(points,data_u,points2,method='cubic').reshape((nx,ny)).T
    interp_v = griddata(points,data_v,points2,method='cubic').reshape((nx,ny)).T

    mag = np.sqrt(interp_v**2 + interp_u**2)
    lw = 0.6#1.2 * mag / np.nanmax(mag) +0.2

    print(np.nanmax(mag))

    ax = axes[i]

    # cnt = ax.pcolormesh(np.degrees(lons),
    #                   np.degrees(lats),
    #                   mag,
    #                 #   levels=levels,
    #                   norm=norm,
    #                   vmin=0,
    #                   vmax=1.)
    cnt = ax.contourf(np.degrees(lons),
                      np.degrees(lats),
                      mag,
                      edgecolor='k',
                      levels=levels,
                      norm=norm,
                      vmin=dmin,
                      vmax=dmax)

    grey = 0.9
    stream  = ax.streamplot(np.degrees(lons),
                             np.degrees(lats),
                             interp_u,
                             interp_v,
                             color=(grey,grey,grey),
                            #  cmap=plt.cm.viridis,
                             norm=norm,
                            #  cmap=cmocean.cm.haline,
                             linewidth=lw,
                             arrowsize=0.8,
                             arrowstyle='fancy',
                             density=1.3)


    # ax.annotate("$t/T = \\num{0.5}$",[340,80])

    ax.annotate("$t/T = \\num{" + "{:.2f}".format(times[i]) + "}$",
                xy=(0.86, 1.0),
                xytext=(0, 0),
                xycoords='axes fraction',
                textcoords='offset points',
                size=10, ha='center', va='bottom')

    ax.xaxis.set_ticks(xticks)
    ax.yaxis.set_ticks(yticks)

    ax.set_xlim(np.degrees([0,2*np.pi]))
    ax.set_ylim(np.degrees([-np.pi*0.5,np.pi*0.5]))

    # ax.set_facecolor((0.95, 0.95, 0.95))



# divider = make_axes_locatable(axes[-1])
# cax = divider.append_axes("right", size="5%", pad=0.8)

fig.subplots_adjust(right=0.8)
cax = fig.add_axes([0.82, 0.25, 0.02, 0.5])
# fig.colorbar(axes[-1], cax=cbar_ax)
ticks = np.linspace(dmin,dmax,5)
clb = plt.colorbar(cnt,
                   cax=cax,
                   orientation='vertical',
                   label="Speed [\si{\metre\per\second}]",
                   norm=norm,
                   ticks=ticks)

for ax in axes:
   ax.set_xlim(np.degrees([0,2*np.pi]))
   ax.set_ylim(np.degrees([-np.pi*0.5,np.pi*0.5]))

# ax
   # ax.set_aspect('equal')
# clb = plt.colorbar(stream.lines,
#                     ax=cax,
#                     fraction=0.06,
#                     pad=0.15,
#                     orientation='horizontal',
#                     label="Velocity [\si{\metre\per\second}]")

# axes[-1].set_title("Enceladus Eccentricity Tide Flow Field")

for ax in l_axes:
    ax.set_ylabel("Latitude [\si{\degree}]")

for ax in b_axes:
    ax.set_xlabel("Longitude [\si{\degree}]")
# axes[3-1].set_xlabel("Longitude [\si{\degree}]")
# axes[4-1].set_xlabel("Longitude [\si{\degree}]")

plt.subplots_adjust(wspace=0.1, hspace=0.15)

# fig.savefig("inf_lid_ecc_enc.pdf")
fig.savefig("inf_gan_ecc_west.png",bbox_inches='tight',dpi=300)


plt.show()

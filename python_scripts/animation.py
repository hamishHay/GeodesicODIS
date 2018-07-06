#!/home/hamish/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
import sys
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpe
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-s","--start",
                    help="Not yet",
                    type=int,
                    default=0)
parser.add_argument("-e","--end",
                    help="Not yet",
                    type=int,
                    default=10)
parser.add_argument("-f","--filename",
                    help="Name of the savefile",
                    type=str,
                    default="animation_odis.mp4")
parser.add_argument("-t","--title",
                    help="Animation title",
                    type=str,
                    default="ODIS Animation")
parser.add_argument("-p","--plot",
                    help="Plot type",
                    type=str,
                    default="velocity")
parser.add_argument("-h0","--h0",
                    help="Ocean thickness",
                    type=float,
                    default=1000.0)
parser.add_argument("-q","--scale",
                    help="velocity quiver scale",
                    type=float,
                    default=1e-2)

args = parser.parse_args()

start = args.start
stop = args.end
filename = args.filename
title = args.title
p_type = args.plot
h = args.h0
scale = args.scale

P = 15.95

dt = P * 0.25
time = 0.0


print(p_type)
# THESE TWO LINES ARE REQD TO PRODUCE SAVEABLE MOVIE
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=4000)

in_file = h5py.File("DATA/data.h5", 'r')

DATA_V = np.array(in_file["north velocity"])[start:stop+1,:,:]
DATA_U = np.array(in_file["east velocity"])[start:stop+1,:-1,:]

DATA_V = DATA_V[:,::-1,:]
DATA_U = DATA_U[:,::-1,:]

DATA = DATA_V

if p_type == "velocity":
    DATA_MAG = np.sqrt(DATA_U**2 + DATA_V**2)
elif p_type == "dissipation":
    DATA_MAG = 1e3 * h * 1e-6 * (DATA_U**2 + DATA_V**2)
elif p_type == 'displacement':
    DATA_MAG = np.array(in_file["displacement"])[start:stop+1,:,:]

frames = stop-start
vmin = np.amin(DATA_MAG)
vmax = np.amax(DATA_MAG)

levels = np.linspace(vmin,vmax,12)

lon = np.linspace(0,360,np.shape(DATA)[2])
lat = np.linspace(-90,90,np.shape(DATA)[1])

LON,LAT = np.meshgrid(lon,lat)


fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

ax.set_xlim([0,360])
ax.set_ylim([90,-90])

if p_type == "velocity":
    qv = ax.quiver(lon[::3],lat[::3],DATA_U[0,::3,::3],DATA_V[0,::3,::3],scale=scale,pivot='mid',color='k')
    ct = ax.imshow(DATA_MAG[0,:,:],vmin=vmin,vmax=vmax,cmap=plt.cm.viridis,interpolation='bicubic',extent=[0,360,90,-90])
elif p_type == "dissipation":
    ct = ax.imshow(DATA_MAG[0,:,:],vmin=vmin,vmax=vmax,cmap=plt.cm.plasma,interpolation='bicubic',extent=[0,360,90,-90])
    ct2 = ax.contour(lon,lat,DATA_MAG[0,:,:],levels=levels,colors='k')
elif p_type == "displacement":
    ct = ax.imshow(DATA_MAG[0,:,:],vmin=vmin,vmax=vmax,cmap=plt.cm.plasma,interpolation='bicubic',extent=[0,360,90,-90])
    ct2 = ax.contour(lon,lat,DATA_MAG[0,:,:],levels=levels,colors='k')

#ax.clear()

ax.set_aspect('equal')
ax.invert_yaxis()
ax.set_ylabel("Latitude (deg)")
ax.set_xlabel("Longitude (deg)")
plt.suptitle(title)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
c = plt.colorbar(ct,cax=cax)

#ct = ax.imshow((sum(DATA_MAG[:,:,:])/len(DATA_MAG)),cmap=plt.cm.plasma,extent=[0,360,90,-90])
#ct2 = ax.contour(lon,lat,sum(DATA_MAG[:,:,:])/len(DATA_MAG))
#plt.show()

if p_type == "velocity":
    c.set_label("Velocity [m/s]")
elif p_type == "dissipation":
    c.set_label("Dissipated Energy Flux [W/m^2]")

def animate(i):
    global time
    if p_type == "velocity":
        ct.set_array(DATA_MAG[i,:,:])
        qv.set_UVC(DATA_U[i,::3,::3],DATA_V[i,::3,::3])
    elif p_type == "dissipation" or p_type == "displacement":
        global ct2
        ct.set_array(DATA_MAG[i,:,:])
        for contour in ct2.collections:
            contour.remove()
        ct2 = ax.contour(lon,lat,DATA_MAG[i,:,:],levels=levels,colors='k')


    time += dt
    ax.set_title('t = %.1f orbits'%(time/P))
    print("Rendering frame ", i)
    return ax

interval = 500#in seconds
ani = animation.FuncAnimation(fig,animate,frames=frames,interval=50,blit=False)

plt.show()

time = 0.0

ani.save(filename,writer=writer)

def get_area(lon,lat):
    dlon = np.deg2rad(lon[1]-lon[0])
    area = np.zeros(len(lat))
    for i in range(len(lat)-1):
        area[i] = (2575e3)**2 * (-np.cos(np.deg2rad(lat[i+1])) + np.cos(np.deg2rad(lat[i])))*dlon
    return area

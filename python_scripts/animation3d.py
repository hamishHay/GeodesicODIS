import numpy as np
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = '/bin/ffmpeg'
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as manimation
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import FieldPlot
import os
import types
from mpl_toolkits.mplot3d import Axes3D
import sys

# plt.rc('font', family='serif')
# plt.rc('font',serif='Palatino')
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rc('text', usetex=True)
# plt.rc('text.latex',preamble='\\usepackage{siunitx}')

plt.rcParams['contour.negative_linestyle'] = 'solid'



# plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
directory = os.getcwd()
file = "/home/hamish/ODIS/test_eta"
# file = "/inResTime"

def GrabData(i,name_id):
    num = i
    snum = str(num)
    name = [file+"/EastVelocity/u_vel_"+snum+".txt", file+"/NorthVelocity/v_vel_"+snum+".txt", file+ "/Displacement/eta_"+snum+".txt"]


    # path = directory + name[name_id]
    path = name[name_id]
    data = []
    if os.path.exists(path):
        print("Retrieving data from:", path)
        init = open(path,'r')
        lines = init.readlines()
        init.close()

        for line in lines:
            line = line.split('\t')
            data.append([float(j) for j in line[:-1]])
        data = np.array(data)

        return data

def sph2cart(r,theta,phi):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z



first = 1
last = int(sys.argv[1])
frames = range(first,last)

out_time = 0.01
res = 10

data_diss = np.loadtxt(file+"/Energy/diss_energy.txt")[first:last]
orbit_num = np.linspace((first-1)*out_time,len(data_diss)*out_time,len(data_diss))

print(orbit_num[-1])

mins = []
maxs = []

for i in frames:
    data = GrabData(i,2)
    mins.append(data.min())
    maxs.append(data.max())

disp_min = min(mins)
disp_max = max(maxs)

mins = []
maxs = []

# for i in frames:
#     data_u = GrabData(i,0)[:-1,:]
#     data_v = GrabData(i,1)
#     data = np.sqrt(data_u**2 + data_v**2)
#     mins.append(data.min())
#     maxs.append(data.max())
#
# vel_min = min(mins)
# vel_max = max(maxs)


data_disp = GrabData(10,2)
# data_u = GrabData(1,0)[:-1,:]
# data_v = GrabData(1,1)



#plt.plot(orbit_num,data_diss)
#plt.show()
#fig = plt.gcf()
#fig.savefig("test.png")

lon = np.linspace(0,360,len(data_disp[0]))
lat = np.linspace(90,-90,len(data_disp))

colat = np.deg2rad(90.0 - lat)
lon = np.deg2rad(lon)
data_disp = data_disp + 250e3

Lon, Colat = np.meshgrid(lon, colat)

print(data_disp.shape, Colat.shape, Lon.shape)
X, Y, Z = sph2cart(data_disp, Colat, Lon)

# x2 = np.linspace(0,360,len(data_v[0]))
# y2 = np.linspace(90,-90,len(data_v))




# X, Y = np.meshgrid(x, y)


cont_num = 21
levels = np.linspace(disp_min,disp_max,cont_num)
# levels_u = np.linspace(vel_min,vel_max,cont_num)
cmap = cm.gist_ncar

stride = 1

fig = plt.figure(figsize=(20/1.2,8/1.2),dpi=300)
# fig = plt.figure()
ax1 = plt.subplot(111, projection='3d',aspect='equal')
quad = ax1.plot_wireframe(X, Y, Z, rstride=stride, cstride=stride)

# plt.show()

# print(0.05*vel_max/vel_max)

X, Y, Z = sph2cart((16*disp_max)**2, Colat, Lon)
mag_max = 16*disp_max#np.sqrt(X**2 + Y**2 + Z**2)

X, Y, Z = sph2cart((disp_min + 15*disp_max)**2, Colat, Lon)
mag_min = disp_min + 15*disp_max#np.sqrt(X**2 + Y**2 + Z**2)


X2 = X**2
Y2 = Y**2
def init():
    quad = ax1.plot_wireframe(X, Y, data_disp, rstride=stride, cstride=stride)
    # quad = ax1.plot_wireframe(X, Y, data_disp, rstride=stride, cstride=stride)

    return quad,


def animate(i):
    # fig.clf()
    ax1.cla()
    data_disp = GrabData(i,2)
    r = data_disp + 15*disp_max
    X, Y, Z = sph2cart(r, Colat, Lon)
    mag = np.sqrt(X**2 + Y**2 + Z**2)
    N = (mag - mag_min)/(mag_max - mag_min)

    # print(N)

    # quad = ax1.plot_wireframe(X, Y, Z, rstride=stride, cstride=stride)
    quad = ax1.plot_surface(X, Y, Z, rstride=stride, cstride=stride,facecolors=cmap(N),shade=False)

    ax1.set_zlim([-18*disp_max,18*disp_max])
    ax1.set_xlim([-18*disp_max,18*disp_max])
    ax1.set_ylim([-18*disp_max,18*disp_max])
    # ax2.plot(seeds[0], seeds[1], 'ro')
    plt.hold('off')


    return quad,

print(frames)

anim = manimation.FuncAnimation(fig, animate,init_func=init,frames = frames,blit=False,interval=1000,repeat=False)
mywriter = manimation.FFMpegWriter(fps=20,bitrate=8000)

# anim.save('h1000_cD0.1_2.gif', writer = 'imagemagick',dpi=300,fps=10)
anim.save('var_thick.mp4', writer = mywriter,dpi=300, extra_args=['-vcodec', 'libx264'])

import numpy as np
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = '/bin/ffmpeg'
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as manimation
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import FieldPlot
import os
import types
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

for i in frames:
    data_u = GrabData(i,0)[:-1,:]
    data_v = GrabData(i,1)
    data = np.sqrt(data_u**2 + data_v**2)
    mins.append(data.min())
    maxs.append(data.max())

vel_min = min(mins)
vel_max = max(maxs)


data_disp = GrabData(1,2)
data_u = GrabData(1,0)[:-1,:]
data_v = GrabData(1,1)



#plt.plot(orbit_num,data_diss)
#plt.show()
#fig = plt.gcf()
#fig.savefig("test.png")

x = np.linspace(0,360,len(data_disp[0]))
y = np.linspace(90,-90,len(data_disp))

x2 = np.linspace(0,360,len(data_v[0]))
y2 = np.linspace(90,-90,len(data_v))

X, Y = np.meshgrid(x, y)


cont_num = 21
levels = np.linspace(disp_min,disp_max,cont_num)
levels_u = np.linspace(vel_min,vel_max,cont_num)
cmap = cm.winter


fig = plt.figure(figsize=(20/1.2,8/1.2),dpi=300)
# fig = plt.figure()
gs = gridspec.GridSpec(2,2,height_ratios=[4,1])
ax1 = plt.subplot(gs[0,0], aspect='equal')
quad = ax1.pcolormesh(x, y, data_disp,cmap=cmap,vmin=disp_min,vmax=disp_max)
plt.hold('on')
quad2 = ax1.contour(x,y,data_disp,levels=levels,colors='k',linewidths=0.4)

# pos1 = ax1.get_position()
# axColor = plt.axes([pos1.x0*1.05 + pos1.width * 1.05, pos1.y0, 0.01, pos1.height])

# cbpos1 = cb1.get_position()
# cb1.set_position

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.06)
cb1 = fig.colorbar(quad,cax=cax1)

ax2 = plt.subplot(gs[0,1], aspect='equal')

quad3 = ax2.pcolormesh(x, y, data_u,cmap=cmap,vmin=vel_min,vmax=vel_max)
plt.hold('on')
# quad5 = ax2.contour(x,y,data_u,levels=levels_u,colors='k',linewidths=0.4)
# quad4 = ax2.streamplot(x2, y2, data_u,data_v,color='k',linewidth=0.4,density=4)#, start_points=seeds)
quad4 = ax2.quiver(x2[::res], y2[::res], data_u[::res,::res],data_v[::res,::res],color='k')#,scale=0.05*np.amax(vel)/vel_max)#,linewidth=0.4,density=0.8,arrowsize=3)#, start_points=seeds)


plt.hold('off')

divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="5%", pad=0.06)
cb2 = fig.colorbar(quad3,cax=cax2)

ax3 = plt.subplot(gs[1,:])
line1, = ax3.plot([],[],'r-',linewidth=1.0)

ax3.set_xlim(0, orbit_num[-1])
ax3.set_ylim(0, 50)#1.2*max(data_diss))
ax3.set_ylabel('Dissipation (\si{\watt\per\metre\squared})')
ax3.set_xlabel('Orbit Number')
#ax3.set_yscale('log')
# plt.show()

print(0.05*vel_max/vel_max)

def init():
    quad = ax1.pcolormesh(x, y, data_disp,cmap=cmap,vmin=disp_min,vmax=disp_max)
    quad2 = ax1.contour(x,y,data_disp,levels=levels,colors='k',linewidths=0.7)

    vel = np.sqrt(data_u**2 + data_v**2)
    quad3 = ax2.pcolormesh(x, y, vel,cmap=cmap,vmin=vel_min,vmax=vel_max)
    plt.hold('on')
    quad4 = ax2.streamplot(x2, y2, data_u,data_v,color='k',linewidth=0.4,density=0.8,arrowsize=3)#, start_points=seeds)
    plt.hold('off')

    line1, = ax3.plot([],[])

    return line1, quad, quad3,


def animate(i):
    # fig.clf()
    ax1.cla()
    data_disp = GrabData(i,2)
    quad = ax1.pcolormesh(x, y, data_disp,cmap=cmap,vmin=disp_min,vmax=disp_max)
    quad2 = ax1.contour(x,y,data_disp,levels=levels,colors='k',linewidths=0.7)
    ax1.set_xlim(x[0], x[-1])
    ax1.set_ylim(y[0], y[-1])
    ax1.set_yticks([90,60,30,0,-30,-60,-90])
    ax1.set_xticks(np.linspace(0,360,9))
    ax1.invert_yaxis()

    ax1.set_ylabel("Latitude ($^{\circ}$)")
    ax1.set_xlabel('Longitude ($^{\circ}$)')
    ax1.set_title("Displacement, $\\eta$ (m)")

    ax2.cla()
    data_u = GrabData(i,0)[:-1,:]
    data_v = GrabData(i,1)
    vel = np.sqrt(data_u**2 + data_v**2)
    quad3 = ax2.pcolormesh(x, y, vel,cmap=cmap,vmin=vel_min,vmax=vel_max)
    quad5 = ax2.contour(x2,y2,vel,levels=levels_u,colors='k',linewidths=0.7)
    plt.hold('on')
    vel = vel/vel_max
    quad4 = ax2.quiver(x2[::res], y2[::res], data_u[::res,::res],data_v[::res,::res],color='k',scale=35)#,linewidth=0.4,density=0.8,arrowsize=3)#, start_points=seeds)
    # quad4.set_UVC(data_u[::res,::res],data_v[::res,::res])


    # ax2.plot(seeds[0], seeds[1], 'ro')
    plt.hold('off')
    ax2.set_xlim(x[0], x[-1])
    ax2.set_ylim(y[0], y[-1])
    ax2.set_yticks([90,60,30,0,-30,-60,-90])
    ax2.set_xticks(np.linspace(0,360,9))
    ax2.set_title("Velocity, $\\left|\\mathbf{u}\\right|$ (\\si{\\metre\\per\\second})")
    ax2.invert_yaxis()


    ax2.set_ylabel("Latitude ($^{\circ}$)")
    ax2.set_xlabel("Longitude ($^{\circ}$)")

    line1.set_data(orbit_num[:i], data_diss[:i])


    return line1, quad, quad3,

gs.tight_layout(fig)
print(frames)

anim = manimation.FuncAnimation(fig, animate,init_func=init,frames = frames,blit=False,interval=1000,repeat=False)
mywriter = manimation.FFMpegWriter(fps=20,bitrate=8000)

# anim.save('h1000_cD0.1_2.gif', writer = 'imagemagick',dpi=300,fps=10)
anim.save('var_thick.mp4', writer = mywriter,dpi=300, extra_args=['-vcodec', 'libx264'])

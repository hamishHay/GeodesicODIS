# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import numpy as np
import math
import h5py
from mpl_toolkits.basemap import Basemap
import sys
import matplotlib.ticker as ticker

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

grid = np.loadtxt("input_files/grid_l6.txt",skiprows=1,usecols=(1,2))

x = grid[:,1]
y = grid[:,0]

n = int(sys.argv[1])

in_file = h5py.File("DATA/data.h5", 'r')
# in_file = h5py.File("DATA_G6_T1000/data.h5", 'r')

data_eta = np.array(in_file["displacement"][-100:-1])

p_spole = data_eta[:,1]
# p_spt_avg = np.mean(data_eta[:,y<-60.0],axis=1)

data_eta = [data_eta[0], data_eta[25], data_eta[50], data_eta[75]]

p_max = np.amax(data_eta)/1e6
p_min = np.amin(data_eta)/1e6

print("TRIANGULATING POSITIONS")

triang = tri.Triangulation(x, y)

fig, axes = plt.subplots(ncols=2, nrows=2, sharey=True, sharex=True, figsize=(6.5,4))
axes = [axes[0][0], axes[0][1], axes[1][0], axes[1][1]]

levels = np.linspace(-0.12,0.12,9)
t = [0.0,0.25,0.5,0.75]

for i in range(4):
    tcnt = axes[i].tricontourf(triang, data_eta[i]/1e6, levels=levels, cmap = plt.cm.coolwarm, vmin=p_min, vmax=p_max)
    cnt = axes[i].tricontour(triang, data_eta[i]/1e6, levels=levels, colors='k', linewidths=0.5)
    # axes[i].set_title("$t/T =" + str(t[i]) + "$")
    axes[i].text(5,-85,"$t/T =" + str(t[i]) + "$", size=8)
    # cbar.set_title("Pressure [\si{\mega\pascal}]")

# cb = plt.colorbar()
# labels = np.arange(0,k,1)
# loc    = labels + .5
# cb.set_ticks(loc)
# cb.set_ticklabels(labels)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cb = fig.colorbar(tcnt, cax=cbar_ax, label='Pressure [\si{\mega\pascal}]')

cb.set_ticks([-0.12,-0.09,-0.06,-0.03,0.0,0.03, 0.06,0.09,0.12])

axes[0].set_ylabel("Latitude [\si{\degree}]")
axes[2].set_ylabel("Latitude [\si{\degree}]")
axes[2].set_xlabel("Longitude [\si{\degree}]")
axes[3].set_xlabel("Longitude [\si{\degree}]")

fig.suptitle("Pressure excess in Enceladus' ocean")

fig.savefig("p_excess.pdf", bbox_inches='tight')

# cbar = plt.colorbar(tcnt, ax = axes[4], orientation='vertical',label="Pressure [\si{\mega\pascal}]")
fig, ax = plt.subplots(figsize=(6,4))

time = np.linspace(0,1,len(p_spole))

ax.plot(time, p_spole/1e6, lw=1.5)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.set_xlabel("Time")
ax.set_ylabel("Pressure excess [\si{\mega\pascal}]")
ax.grid(which='both')
ax.set_title("Pressure excess at Enceladus' south pole")
ax.set_ylim([0.0390, 0.0420])

fig.savefig("p_excess_spt.pdf", bbox_inches='tight')

plt.show()

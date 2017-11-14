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
x = np.radians(grid[:,1])
y = np.radians(grid[:,0])
triang = tri.Triangulation(x, y)

grad = False
div = True

if div == True:
    data = np.loadtxt("error_n1.txt").T

    anlyt = data[4,:]
    num = data[5]



    err = np.sqrt((anlyt - num)**2.0)
    err2 = abs(anlyt-num)/abs(anlyt)*100.0

    for i in range(len(x)):
        print(anlyt[i],'\t\t',num[i],'\t\t',err2[i])

    print(np.nanmean(err2))

    # err = abs(anlyt-num)/abs((anlyt+num)/2)*100 #abs(anlyt-num)/anlyt * 100.0

    # err = np.log10(err)

    cmax = max(np.amax(anlyt), np.amax(num))
    cmin = min(np.amin(anlyt), np.amin(num))

    levels = np.linspace(cmin,cmax,11)

    print("TRIANGULATING POSITIONS")
    # Create the Triangulation; no triangles so Delaunay triangulation created.

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(12,3),dpi=120)

    tcnt = ax1.tricontourf(triang, anlyt, levels=levels, cmap = plt.cm.winter)
    cnt = ax1.tricontour(triang, anlyt, levels=levels, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
    ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \cdot \\vec{u}$")

    tcnt = ax2.tricontourf(triang, num, levels=levels, cmap = plt.cm.winter)
    cnt = ax2.tricontour(triang, num, levels=levels, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
    ax2.set_title("Numerical, $s_{num} = \\nabla \cdot \\vec{u}$")

    tcnt = ax3.tricontourf(triang, err, 11, cmap = plt.cm.winter)
    cnt = ax3.tricontour(triang, err, 11, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
    ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

    fig.suptitle("Divergence Operator Test, $n=1, m=1$", y=1.05)

    fig.savefig("div_test_n1_m1.pdf", bbox_inches='tight')
# plt.tight_layout()

if grad==True:
    data = np.loadtxt("error_grad_n1.txt").T

    anlyt1 = data[0,:]
    anlyt2 = data[1,:]

    num1 = data[2,:]
    num2 = data[3,:]

    err1 = np.sqrt((anlyt1 - num1)**2.0)
    err2 = np.sqrt((anlyt2 - num2)**2.0)

    err3 = abs(anlyt1-num1)/abs(anlyt1)*100.0
    err4 = abs(anlyt2-num2)/abs(anlyt2)*100.0

    for i in range(len(x)):
        print(anlyt2[i],'\t\t',num2[i],'\t\t',err4[i])

    err4[np.isinf(err4)] = np.nan
    print(np.nanmean(err4))

    cmax1 = max(np.amax(anlyt1), np.amax(num1))
    cmin1 = min(np.amin(anlyt1), np.amin(num1))

    levels1 = np.linspace(cmin1,cmax1,11)

    cmax2 = max(np.amax(anlyt2), np.amax(num2))
    cmin2 = min(np.amin(anlyt2), np.amin(num2))

    print(cmin1, cmax2)

    levels2 = np.linspace(cmin2,cmax2,11)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2,ncols=3,figsize=(12,6),dpi=120)

    tcnt = ax1.tricontourf(triang, anlyt1, levels=levels1, cmap = plt.cm.winter)
    cnt = ax1.tricontour(triang, anlyt1, levels=levels1, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
    ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \\beta$")

    tcnt = ax2.tricontourf(triang, num1, levels=levels1, cmap = plt.cm.winter)
    cnt = ax2.tricontour(triang, num1, levels=levels1, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
    ax2.set_title("Numerical, $s_{num} = \\nabla \\beta$")

    tcnt = ax3.tricontourf(triang, err1, 11, cmap = plt.cm.winter)
    cnt = ax3.tricontour(triang, err1, 11, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
    ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

    tcnt = ax4.tricontourf(triang, anlyt2, levels=levels2, cmap = plt.cm.winter)
    cnt = ax4.tricontour(triang, anlyt2, levels=levels2, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax4, orientation='horizontal')
    ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \\beta$")

    tcnt = ax5.tricontourf(triang, num2, levels=levels2, cmap = plt.cm.winter)
    cnt = ax5.tricontour(triang, num2, levels=levels2, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax5, orientation='horizontal')
    ax2.set_title("Numerical, $s_{num} = \\nabla \\beta$")

    tcnt = ax6.tricontourf(triang, err2, 11, cmap = plt.cm.winter)
    cnt = ax6.tricontour(triang, err2, 11, colors='k', linewidths=0.5)
    cbar = plt.colorbar(tcnt, ax = ax6, orientation='horizontal')
    ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

    fig.suptitle("Gradient Operator Test, $n=2, m=2$")

# fig.savefig("grad_test_n1_m1.pdf", bbox_inches='tight')


plt.show()

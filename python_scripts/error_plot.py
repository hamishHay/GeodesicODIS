# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import numpy as np
import math
import h5py
import sys
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import griddata

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

fname = "div_test.txt"

infile = np.loadtxt("div_test.txt")
print(infile)

grid = np.loadtxt("input_files/grid_l7.txt",skiprows=1,usecols=(1,2))
x = np.degrees(grid[2:,1])
y = np.degrees(grid[2:,0])
triang = tri.Triangulation(x, y)

grad = True
div = False

yi = np.arange(-90,90,4.0)
xi = np.arange(0,360,4.0)
xi,yi = np.meshgrid(xi,yi)

# set mask
# mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)

# interpolate

# data = np.loadtxt(fname,skiprows=21).T

anlyt1 = infile[2:,1]
num1 = infile[2:,2]

print(len(x), len(num1))

# zi = griddata((x,y),num1,(xi,yi),method='linear')
# err = np.sqrt((anlyt - num)**2.0)
# err = np.sqrt((anlyt - num)**2.0)
#
# err2 = abs(anlyt-num)/abs(anlyt)*100
#
# err2 = np.log10(err2)
# err2[np.isinf(err2)] = np.nan
# err5 = err2[~np.isnan(err2)]
# # err = (anlyt - num)**2.0
#
# # plt.hist(err2,bins=30)
# # ax = plt.gca()
# # ax.set_yscale('log')
#
#
# cmax = max(np.amax(anlyt), np.amax(num))
# cmin = min(np.amin(anlyt), np.amin(num))
#
# levels = np.linspace(cmin,cmax,11)
# # err_levels = np.linspace(-4, 5,10)
#
# print("TRIANGULATING POSITIONS")
# # Create the Triangulation; no triangles so Delaunay triangulation created.
#
# fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(12,3),dpi=120)
#
# tcnt = ax1.tricontourf(triang, anlyt, levels=levels, cmap = plt.cm.winter)
# cnt = ax1.tricontour(triang, anlyt, levels=levels, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
# ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \cdot \\vec{u}$")
#
# tcnt = ax2.tricontourf(triang, num, levels=levels, cmap = plt.cm.winter)
# cnt = ax2.tricontour(triang, num, levels=levels, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
# ax2.set_title("Numerical, $s_{num} = \\nabla \cdot \\vec{u}$")
#
# tcnt = ax3.tricontourf(triang, err, 9, cmap = plt.cm.winter)
# cnt = ax3.tricontour(triang, err, 9, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
# ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")
#
# fig.suptitle("Divergence Operator Test, $n=1, m=1$", y=1.05)
#
# # fig.savefig("div_test_n1_m1.pdf", bbox_inches='tight')
# # plt.tight_layout()
#
#
# data = np.loadtxt(fname,skiprows=21).T
#
# anlyt1 = data[0,2:]
# anlyt2 = data[1,2:]
#
# num1 = data[2,2:]
# num2 = data[3,2:]

# err1 = (anlyt1 - num1)**2.0
# err2 = (anlyt2 - num2)**2.0


# err1 = np.sqrt((anlyt1 - num1)**2.0)
# err2 = np.sqrt((anlyt2 - num2)**2.0)

# err1 = (anlyt1 - num1)
# # err2 = (anlyt2 - num2)
#
# err3 = abs(anlyt1-num1)/abs(anlyt1)*100
# # err4 = abs(anlyt2-num2)/abs(anlyt2)*100
#
# err3 = np.log10(err3)
# err3[np.isinf(err3)] = np.nan
# err3 = err3[~np.isnan(err3)]

# print(err3)
#
# err4 = np.log10(err4)
# err4[np.isinf(err4)] = np.nan
# err4 = err4[~np.isnan(err4)]

# fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(10,3),dpi=120)
#
#
# (ax1, ax2, ax3) = axes
# ax1.hist(err5,bins=20,normed=1,edgecolor='k',lw=0.5)
# ax2.hist(err3,bins=20,normed=1,edgecolor='k',lw=0.5)
# ax3.hist(err4,bins=20,normed=1,edgecolor='k',lw=0.5)
#
# ax1.set_title("Divergence Operator")
# ax2.set_title("Gradient Operator (x)")
# ax3.set_title("Gradient Operator (y)")
#
# ml = MultipleLocator(1)
#
# for ax in axes:
#     ax.set_xlim([-4,4])
#     ax.set_xlabel("$\log_{10} E$ [\si{\percent}]")
#     ax.xaxis.set_minor_locator(ml)
# # ax = plt.gca()
# # fig.savefig("error_dist2.pdf",bbox_inches='tight')
#
# cmax1 = max(np.amax(anlyt1), np.amax(num1))
# cmin1 = min(np.amin(anlyt1), np.amin(num1))
#
# levels1 = np.linspace(cmin1,cmax1,11)
#
# cmax2 = max(np.amax(anlyt2), np.amax(num2))
# cmin2 = min(np.amin(anlyt2), np.amin(num2))
#
# print(cmin1, cmax2)
#
# levels2 = np.linspace(cmin2,cmax2,11)

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2,ncols=3,figsize=(12,6),dpi=120)


# ax1.scatter(x, y)
# tcnt = ax1.contourf(xi, yi, zi, cmap = plt.cm.winter)
# cnt = ax1.contour(xi, yi, zi, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
# ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \\beta$")
#
# tcnt = ax2.contourf(xi, yi, zi, cmap = plt.cm.winter)
# cnt = ax2.contour(xi, yi, zi, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
# ax2.set_title("Numerical, $s_{num} = \\nabla \\beta$")
#
# tcnt = ax3.contourf(xi, yi, zi, 9, cmap = plt.cm.winter)
# cnt = ax3.contour(xi, yi, zi, 9, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
# ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

tcnt = ax1.tricontourf(triang, anlyt1, cmap = plt.cm.winter)
cnt = ax1.tricontour(triang, anlyt1, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax1, orientation='horizontal')
ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \\beta$")

tcnt = ax2.tricontourf(triang, num1, cmap = plt.cm.winter)
cnt = ax2.tricontour(triang, num1, colors='k', linewidths=0.5)
cbar = plt.colorbar(tcnt, ax = ax2, orientation='horizontal')
ax2.set_title("Numerical, $s_{num} = \\nabla \\beta$")

tcnt = ax3.tricontourf(triang, abs(anlyt1-num1), 9, cmap = plt.cm.winter)
cnt = ax3.tricontour(triang, abs(anlyt1-num1), 9, colors='k', linewidths=0.1)
cbar = plt.colorbar(tcnt, ax = ax3, orientation='horizontal')
ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

# ax1.triplot(triang)

# tcnt = ax4.tricontourf(triang, anlyt2, levels=levels2, cmap = plt.cm.winter)
# cnt = ax4.tricontour(triang, anlyt2, levels=levels2, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax4, orientation='horizontal')
# ax1.set_title("Analytical,  $s_{anlyt} = \\nabla \\beta$")
#
# tcnt = ax5.tricontourf(triang, num2, levels=levels2, cmap = plt.cm.winter)
# cnt = ax5.tricontour(triang, num2, levels=levels2, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax5, orientation='horizontal')
# ax2.set_title("Numerical, $s_{num} = \\nabla \\beta$")
#
# tcnt = ax6.tricontourf(triang, err2, 9, cmap = plt.cm.winter)
# cnt = ax6.tricontour(triang, err2, 9, colors='k', linewidths=0.5)
# cbar = plt.colorbar(tcnt, ax = ax6, orientation='horizontal')
# ax3.set_title(" $E = \sqrt{(s_{anylt} - s_{num})^2}$")

# for i in range(len(err2)):
#     ax6.text(x[i], y[i], str(i+2))

# ax6.set_xlim([np.radians(30), np.radians(60)])
# ax6.set_ylim([np.radians(30), np.radians(40)])

fig.suptitle("Gradient Operator Test, $n=1, m=1$")

fig.savefig("/home/hamish/Dropbox/Tests/grad_test_n1_m1.pdf", bbox_inches='tight')


plt.show()

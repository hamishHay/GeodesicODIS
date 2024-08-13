import numpy as np
import h5py
import os
import sys
# from scipy.integrate import simps

import matplotlib as mpl
if "SSH_CONNECTION" in os.environ:
    mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1 import ImageGrid

class mpl_defaults:
    def __init__(self):
        print("PLOTMUN")
        # plt.rc('font', family='serif')
        # plt.rc('font',serif='Palatino')
        # for Palatino and other serif fonts use:
        # rc('font',**{'family':'serif','serif':['Palatino']})
        # plt.rc('text', usetex=True)
        # plt.rc('text.latex',preamble='\\usepackage{siunitx}')

class FieldPlot2:
    def __init__(self, ncols=1, nrows=1):
        self.set_defaults(poster=True)

        self.ncols = ncols
        self.nrows = nrows

        self.bot_ind = nrows - 1
        self.left_ind = 0

        fig = plt.figure(figsize=(ncols*(3*2), nrows*(2*2)))



        self.fig = fig
        self.axes = []

        self.load_simulation_info()
        self.load_grid()

    def plot_scalar(self, data_name='dissipated energy',  slices=[0], b_slice=0, ax_ind=[[0,0]], sample=0, crange=None, lvl_num=11):
        self.im_grid = ImageGrid(self.fig,
                                 111,
                                 nrows_ncols=(self.nrows,self.ncols),
                                 share_all=True,
                                 axes_pad=0.15,
                                 cbar_mode='single',
                                 cbar_location='right')


        if slices==None:
            data = self.load_data(data_name=data_name, avg=True, avg_lim=[-100*100,-1])
            slices = [0]
            data = [data[2:]*1000]
        else:
            data = self.load_data(data_name=data_name, slices=slices, base_slice=b_slice)

        if crange == None:
            cmin = np.amin(data)
            cmax = np.amax(data)
        else:
            cmin = crange[0]
            cmax = crange[1]

        clevels = np.linspace(cmin, cmax, lvl_num)
        norm = mpl.colors.Normalize(vmin=cmin,vmax=cmax)

        x = self.grid[:,1]
        y = self.grid[:,0]

        x = x[2:]
        y = y[2:]

        triang = tri.Triangulation(x, y)

        if len(ax_ind) != len(slices):
            raise ValueError("Not enough axes to plot time slices!")

        for i in range(len(slices)):
            ax_current = self.im_grid[i]
            self.axes.append(ax_current)
            cnt = ax_current.tricontourf(triang, data[i], levels=clevels, cmap = plt.cm.plasma)
            ax_current.tricontour(triang, data[i], levels=clevels, colors='k', linewidths=1.0)


        cb = self.im_grid.cbar_axes[0].colorbar(cnt, norm=norm)
        cb.set_label_text(self.data_label[data_name])

        self.set_spatial_plot_defaults(slices,ax_ind)

    def plot_vector(self,  slices=[0], b_slice=0, ax_ind=[[0,0]], sample=0, crange=None, lvl_num=11):
        self.im_grid = ImageGrid(self.fig,
                                 111,
                                 nrows_ncols=(self.nrows,self.ncols),
                                 share_all=True,
                                 axes_pad=0.15,
                                 cbar_mode='single',
                                 cbar_location='right')


        data_name = 'east velocity'
        data_u = self.load_data(data_name=data_name, slices=slices, base_slice=b_slice)

        data_name = 'north velocity'
        data_v = self.load_data(data_name=data_name, slices=slices, base_slice=b_slice)

        data_name = 'displacement'
        data_eta = self.load_data(data_name=data_name, slices=slices, base_slice=b_slice)

        data_name = 'face longitude'
        data_lon = self.load_data(data_name=data_name, slices=None, base_slice=b_slice)

        data_name = 'face latitude'
        data_lat = self.load_data(data_name=data_name, slices=None, base_slice=b_slice)

        data = np.sqrt(data_u**2 + data_v**2)

        print(len(data_lat))

        if crange == None:
            cmin = np.amin(data)
            cmax = np.amax(data)
        else:
            cmin = crange[0]
            cmax = crange[1]

        clevels = np.linspace(cmin, cmax, lvl_num)
        norm = mpl.colors.Normalize(vmin=cmin,vmax=cmax)

        if sample < 0:
            grid2 = self.load_grid(N=self.gg_lvl+sample, return_val=True)
            x = grid2[:,1]
            y = grid2[:,0]

        else:
            x = self.grid[:,1]
            y = self.grid[:,0]

        x1 = data_lon 
        y1 = data_lat

        triang = tri.Triangulation(x, y)
        triang1 = tri.Triangulation(x1, y1)

        # if len(ax_ind) != len(slices):
        #     raise ValueError("Not enough axes to plot time slices!")

        # for i in range(len(slices)):
        #     ax_current = self.im_grid[i]
        #     self.axes.append(ax_current)

        #     # cnt = ax_current.tricontourf(triang1, data[i][:len(x1)], levels=clevels, cmap = plt.cm.plasma)
        #     cnt = ax_current.tricontourf(triang, data_eta[i], levels=np.arange(-0.1, 1.11, 0.1))
        #     print(data_eta[i])
        #     # ax_current.quiver(x1,y1,data_u[:len(x1)],data_v[:len(x1)])#,scale=0.01)

        ax_current = self.im_grid[0]
        self.axes.append(ax_current)

        # cnt = ax_current.tricontourf(triang1, data[i][:len(x1)], levels=clevels, cmap = plt.cm.plasma)
        # cnt = ax_current.tricontourf(triang1, data_u[-1], 11)
        cnt = ax_current.tricontourf(triang1, data_v[-1], 11)
        cnt2 = ax_current.tricontour(triang, data_eta[-1]/1e3, colors='k', linewidths=0.5)#, levels=np.arange(0, 3000, 120))
        # ax_current.tricontour(triang, data_eta[-1], colors='k', linewidths=0.5)#, levels=np.arange(0, 3000, 120))#, levels=np.arange(-0.1, 1.11, 0.1))
        ax_current.clabel(cnt2, cnt2.levels, inline=True, fontsize=6)

        cb = self.im_grid.cbar_axes[0].colorbar(cnt, norm=norm, label="Velocity [m/s]")
        # cb.set_label_text(self.data_label["velocity"])

    

        # self.set_spatial_plot_defaults(slices,ax_ind)

    def plot_avg_vs_time(self,data=None,data_name="dissipation avg output",lims=[0,-1], title=''):

        if data == None:
            data = self.load_data(data_name=data_name, slices=None)
            time = self.load_data(data_name="kinetic avg output", slices=None)/self.period
            vel_u = self.load_data(data_name='east velocity', slices=None)
            vel_v = self.load_data(data_name='north velocity', slices=None)
            eta = self.load_data(data_name='displacement', slices=None)
        # data *= 4*np.pi * (self.radius-)
        P = 2*np.pi/self.omega
        # print(time[-1002], time[-2])
        # time *= P
        data *= 4. * np.pi * (self.radius - (self.h_shell - self.h_ocean))**2.0 /1e9
        # data *= 4. * np.pi * (self.radius)**2.0/1e9
        # data *= (252.1e3 - self.h_shell + self.h_ocean)**2.0
        # data *= (1./(252.1e3))**2.0
        # data *= ((252.1e3-self.h_ocean)/(252.1e3))**3.
        #data *= 4. * np.pi * (self.radius)**2.0/1e9 / (2. * self.drag_coeff)
        # print(data)
        # print(np.mean(data[-102:-2]))
        # print(simps(data[-102:-2], time[-102:-2])/P)
        # print(time[-102]/P, time[-2]/P)
        # print(time[998:1001])

        # Ek = np.sqrt(vel_u**2.0 + vel_v**2.0)

        ax_current = self.fig.add_subplot(111)

        # time = np.arange(self.out_inc,self.t_end+self.out_inc,self.out_inc)
        # # time = time[1:]
        # dt = 20
        # # time = time[-16*dt:]
        # # data = data[-16*dt:]
        # if len(time) != len(data):
        #     data = data[:-1]


        # time = time[-70000*dt:-69800*dt]
        # data = data[-70000*dt:-69800*dt]
        # t2 = np.linspace(0,dt*2,len(data))
        # ax_current.semilogy(data, lw=2.0)

        # print(Ek[-201:-1, 137]/np.amax(Ek[-201:, 137]))


        # print(np.amax(Ek[-102:-1, 137]), np.amin(Ek[-102:-1, 137]))

        # print(np.mean(Ek[-102:-1, 137]))

        # print(time[-1001], time[-1])
        # ax_current.plot(time[-1002:-1], Ek[-1002:-1, 0])

        # ax_current.plot(time[-1002:-1], data[-1002:-1])
        ax_current.semilogy(time[:-1], data[:-1])
        # ax_current


        # ax_current.set_xticks(np.arange(0, t2[-1]+1, 0.5))
        # ax_current.set_xticks(np.arange(0, t2[-1], 0.25), minor=True)

        # ax_current.set_xlim(0,max(t2))
        # ax_current.set_ylim(0.8*np.amin(data[-100:]),1.2*np.amax(data[-100:]))

        # ax_current.set_xlim([0, 3000])
        # ax_current.set_ylim([3e-5, 0.05])

        ax_current.grid(which='major',alpha=0.8)
        ax_current.grid(which='minor',alpha=0.4)

        ax_current.set_xlabel('Time [Orbit \#]')
        ax_current.set_ylabel('Dissipated Power [\si{\giga\watt}]')

        # ax_current.set_xlim([9900, 10000])

        # if title=='':
        #     ax_current.set_title(self.data_title[self.potential])
        # else:
        #     ax_current.set_title(title)
        # from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        # from mpl_toolkits.axes_grid1.inset_locator import mark_inset
        #
        # inset_axes = inset_axes(ax_current,
        #             width="30%", # width = 30% of parent_bbox
        #             height=1., # height : 1 inch
        #             loc=1)

        # print(time[-101:])
        #
        # inset_axes.plot(time[-100:], data[-100:],lw=1.5)
        #
        # inset_axes.grid(which='major',alpha=0.8)
        # inset_axes.grid(which='minor',alpha=0.4)
        #
        # plt.ticklabel_format(style='plain', axis='x')

        plt.show()

    def set_spatial_plot_defaults(self, slices, ax_ind):
        self.im_grid.axes_llc.set_xticks(np.arange(0,350.,45.))
        self.im_grid.axes_llc.set_yticks(np.arange(-90,91.,45.))

        for i in range(len(slices)):
            ax_current = self.im_grid[i]

            ax_current.set_aspect('equal')
            ax_current.set_adjustable('box-forced')
            ax_current.set_ylim([-90,90])
            ax_current.set_xlim([0,360])

            if ax_ind[i][0] == self.bot_ind:
                ax_current.set_xlabel("Longitude [\si{\degree}]")

            if ax_ind[i][1] == self.left_ind:
                ax_current.set_ylabel("Latitude [\si{\degree}]")

    def get_fig(self):
        return self.fig, self.axes

    def load_grid(self,N=0,return_val=False):
        if return_val:
            return np.loadtxt("input_files/grid_l" + str(N) + ".txt",skiprows=1,usecols=(1,2))

        self.grid = np.loadtxt("input_files/grid_l" + str(self.gg_lvl) + ".txt",skiprows=1,usecols=(1,2))

    def load_data(self, data_name='dissipated energy', slices=[0], base_slice=0, avg=False, avg_lim=[0,-1],avg_axis=0, root=''):
        in_file = h5py.File(root + "DATA/data.h5", 'r')

        if slices==None:
            data = in_file[data_name][:]
            return data

        if not avg:
            if base_slice < 0:
                base_slice = len(in_file[data_name]) + base_slice

            slices = [x + base_slice for x in slices]

            data = np.array(in_file[data_name][slices])
        else:
            data = np.mean(np.array(in_file[data_name][avg_lim[0]:avg_lim[1]]),axis=avg_axis)



        return data

    def save_fig(self,name="FieldPlot2_plot.pdf"):
        self.fig.savefig(name, dpi=400, bbox_inches='tight', transparent=False)

    def set_defaults(self, poster=False):
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


        if poster:
            # plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'red', 'ytick.color':'green', 'figure.facecolor':'white'})
            # plt.rc('axes', edgecolor='white')
            # plt.rc('xtick', color='white')
            # plt.rc('ytick', color='white')
            plt.rc('lines', linewidth=2.0)
            #
            plt.rc('xtick', labelsize=12)
            plt.rc('ytick', labelsize=12)



            plt.rc('axes', titlesize=16)
            plt.rc('axes', labelsize=14)
            plt.rc('axes', labelpad=10.0)
            plt.rc('axes', linewidth=1.20)
            print(plt.style.available)
            #plt.style.use(['dark_background'])#, 'presentation'])

# ytick.labelsize : 16

        # plt.rcParams['font.family'] = 'serif'
        # plt.rcParams['font.serif'] = 'Ubuntu'
        # plt.rcParams['font.monospace'] = 'Ubuntu Mono'
        # plt.rc('font', family='serif')
        else:
            plt.rc('lines', linewidth=0.6)

        plt.rc('figure', dpi=120)

        self.data_label = {"dissipated energy": "Dissipated Energy [\si{\milli\watt}]",
                           "velocity": "Speed [\si{\metre\per\second}]",
                           "displacement": "Pressure"}

        self.data_title = {"ECC_EAST": "Eccentricity tide; eastward",
                           "ECC_WEST": "Eccentricity tide; westward",
                           "ECC_RAD":  "Eccentricity tide; radial",
                           "ECC":      "Eccentricity tide",
                           "OBLIQ_EAST": "Obliquity tide; eastward",
                           "OBLIQ_WEST": "Obliquity tide; westward",
                           "OBLIQ":      "Obliquity tide",
                           "FULL":       "Oblquity and Eccentricity tides"}



    def load_simulation_info(self, file_name='input.in'):
        in_file = open(file_name, 'r')

        lines = in_file.readlines()

        for line in lines:
            line = line.split(';')[:-2]
            line[1] = line[1].strip()

            print(line)

            var_name = line[0]
            var_val = line[1]

            if var_name == 'ocean thickness':
                self.h_ocean = float(var_val)

            elif var_name == 'radius':
                self.radius = float(var_val)

            elif var_name == 'friction coefficient':
                self.drag_coeff = float(var_val)

            elif var_name == 'friction type':
                if var_val == 'LINEAR':
                    self.drag_type = "linear drag"
                elif var_val == 'QUADRATIC':
                    self.drag_type = "bottom drag"

            elif var_name == 'potential':
                self.potential = var_val

            elif var_name == 'geodesic grid level':
                self.gg_lvl = int(var_val)

            elif var_name == 'simulation end time':
                self.t_end = float(var_val)

            elif var_name == 'shell thickness':
                self.h_shell = float(var_val)

            elif var_name == 'output time':
                self.out_inc = float(var_val)

            elif var_name == 'angular velocity':
                self.omega = float(var_val)

        self.period = 2.*np.pi/self.omega



if __name__=='__main__':
    # fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4)
    # axes = [[ax1, ax2, ax3, ax4]]
    # ncols = 4
    # nrows = 1
    # Plotter = FieldPlot2(ncols,nrows)
    #
    # plot_ax = [[0,0],[0,1],[0,2],[0,3]]#,
    #         #    [1,0],[1,1],[1,2],[1,3]]
    #
    # slices = [0,25,50,75]
    # # slices = [0,12,25,37,50,62,75,88]
    #
    # Plotter.plot_scalar(data_name='dissipated energy',
    #                     slices=slices,
    #                     b_slice=-100,
    #                     ax_ind=plot_ax,
    #                     crange=None,
    #                     lvl_num=11)

    ncols = 1
    nrows = 1
    Plotter = FieldPlot2(ncols,nrows)

    plot_ax = [[0,0]]#,
            #    [1,0],[1,1],[1,2],[1,3]]

    slices = None
    # slices = [0,12,25,37,50,62,75,88]

    # Plotter.plot_scalar(data_name='dissipated energy',
    #                     slices=slices,
    #                     b_slice=-100,
    #                     ax_ind=plot_ax,
    #                     crange=None,
    #                     lvl_num=11)

    # plt.show()
    iter = int(sys.argv[1])
    Plotter.plot_vector(slices=[iter])#slices=slices,
                        # b_slice=-100,
                        # ax_ind=plot_ax,
                        # crange=None,
                        # lvl_num=11,
                        # sample=-1)

    # slices = np.arange(0,300,1,dtype=np.int)

    # ncols = 1
    # nrows = 1
    # Plotter = FieldPlot2(ncols,nrows)

    # data_d1_u = Plotter.load_data(base_slice=70000,slices=slices,data_name="east velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_lib_east/')
    # data_d2_u = Plotter.load_data(base_slice=70000,slices=slices,data_name="east velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_lib_west/')
    # data_d3_u = Plotter.load_data(base_slice=70000,slices=slices,data_name="east velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_rad/')
    #
    # data_d1_v = Plotter.load_data(base_slice=70000,slices=slices,data_name="north velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_lib_east/')
    # data_d2_v = Plotter.load_data(base_slice=70000,slices=slices,data_name="north velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_lib_west/')
    # data_d3_v = Plotter.load_data(base_slice=70000,slices=slices,data_name="north velocity",root='/home/hamish/Research/InfShell/Ganymede/tests/ecc_rad/')

    # data_u = data_d1_u + data_d2_u# - data_d3_u
    # data_v = data_d1_v + data_d2_v# + data_d3_v




    # data = np.mean(data_u**2 + data_v**2,axis=1)

    # print(data)
    # plt.plot(data)
    # plt.show()
    # print(data)
    # data = np.mean(data, axis=1)
    # Plotter.plot_avg_vs 

    Plotter.save_fig(name='velocity_soln.png')

    fig,ax = plt.subplots()
    data = Plotter.load_data(data_name="dissipation avg output", slices=None)
    print(data)
    print(len(data))
    ax.plot(data[data>0])
    fig.savefig("dissipation.png",dpi=400,bbox_inches='tight')

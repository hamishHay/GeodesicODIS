import numpy as np
import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy as sc
from scipy.interpolate import interp2d
import os
import sys
import h5py

class ODISPlot:
    def __init__(self,orbitnum=0,figdim=(1,2),figsize=(10,3),saveName="ODISPLOT"):
        # Set matplotlib text options
        plt.rc('font', family='serif')
        plt.rc('font',serif='Palatino')
        # for Palatino and other serif fonts use:
        # rc('font',**{'family':'serif','serif':['Palatino']})
        #plt.rc('text', usetex=True)
        #plt.rc('text.latex',preamble='\\usepackage{siunitx}')

        # Set matplotlib rc options
        plt.rcParams['axes.linewidth'] = 0.8

        # Set colourmap
        self.cmap=cm.winter

        # Set label and title sizes
        self.labelSize = 11
        self.slabelSize = 10
        self.titleSize = 12
        self.supTitleSize = 13

        # Set plotting orbit
        self.orbit = orbitnum

        # Set figure dimensions
        self.currentFig = 0
        self.figSize = figsize
        self.figDim = figdim

        # Set figure save name
        self.saveName = saveName

        # Define figure
        self.fig = plt.figure(1,figsize=self.figSize,dpi=120)

    def loadFieldAxisOptions(self,ax):
        plt.rcParams['contour.negative_linestyle'] = 'solid'

        self.xticks=[0,45,90,135,180,225,270,315,360]
        self.yticks=[90,60,30,0,-30,-60,-90]

        ax.xaxis.set_ticks( self.xticks )
        ax.yaxis.set_ticks( self.yticks )

        ax.set_xlabel('Longitude ($^{\circ}$)',fontsize=self.slabelSize)
        ax.set_ylabel('Latitude ($^{\circ}$)',fontsize=self.slabelSize)

        ax.tick_params(axis='both', which='major', labelsize=self.slabelSize)

        ax.axis([0, 360, -90, 90])
        ax.set_aspect('equal')

    def loadFieldColourOptions(self,cb,ax,cticks=[]):
        cb.ax.tick_params(axis='both', which='major', labelsize=self.labelSize)
        cb.outline.set_linewidth(0.7)
        cb.set_label(self.caxisLabel,fontsize=self.labelSize)
        if len(cticks)>0:
            cb.set_ticks(cticks)

    def ReadFieldData(self,field="displacement",directory = os.getcwd()):
        orbitnum=str(self.orbit)

        name = ["/EastVelocity/u_vel_"+orbitnum+".txt", "/NorthVelocity/v_vel_"+orbitnum+".txt", "/Displacement/eta_"+orbitnum+".txt", "/Energy/diss_"+orbitnum+".txt"]

        if field=="displacement":
            name = name[2]
            self.caxisLabel = "Dispclacement, $\\eta$ (m)"
        elif field=="north_vel":
            name = name[1]
            self.caxisLabel = "North Velocity, $v$ (\\si{\\metre\\per\\second})"
        elif field=="east_vel":
            name = name[0]
            self.caxisLabel = "East Velocity, $v$ (\\si{\\metre\\per\\second})"
        elif field=="dissipation":
            name = name[3]
            self.caxisLabel = "Dissipated Energy, [\si{\watt}]"
        else:
            print("Field " + field + "not recognised. Exiting.")
            sys.exit()

        path = directory + name
        self.data = []
        if os.path.exists(path):
            print("Retrieving data from:", path)
            # init = open(path,'r')
            # lines = init.readlines()
            # init.close()
            #
            # for line in lines:
            #     line = line.split('\t')
            #     self.data.append([float(j) for j in line[:-1]])
            # self.data = np.array(self.data)
            self.data = np.loadtxt(path)

        self.x = np.linspace(0,360,len(self.data[0]))
        self.y = np.linspace(90,-90,len(self.data))

        return self.data

    def ReadDissipationData(self,directory = os.getcwd()):
        in_file = h5py.File("DATA/data.h5", 'r')


        self.diss = np.loadtxt(directory) # In Watts
        # self.diss *= 1e3 # Convert to mWatts
        self.resid = abs(self.diss[1:]-self.diss[:-1])*1e-3
        return self.diss

    def PlotField(self,data=[],cformat=2,cticks=[]):
        self.currentFig += 1

        if len(data) == 0:
            data=self.data

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        self.loadFieldAxisOptions(ax)

        p1 = ax.pcolormesh(self.x,self.y,data,cmap=self.cmap,rasterized=False)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        c1 = plt.colorbar(p1, cax = cax, format='%.' +str(cformat) + 'f',orientation="horizontal")

        cs1 = ax.contour(self.x,self.y,data,8,colors='k',linewidths=0.4)

        self.loadFieldColourOptions(c1,ax,cticks)


    def PlotVelocity(self,u,v,cformat=3,scale=1/0.02,cticks=[]):
        self.currentFig += 1
        self.caxisLabel = "Velocity, $\\left|\\mathbf{u}\\right|$ (\\si{\\metre\\per\\second})"

        u=u[:-1][:]

        mag = np.sqrt(u**2+v**2)

        x = np.linspace(0,360,len(u[0]))
        y = np.linspace(-90,90,len(u))

        mag = np.flipud(mag)

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        p1 = ax.pcolormesh(x,y,mag,cmap=self.cmap,vmin=0,rasterized=False)
        cs1 = ax.contour(x,y,mag,8,colors='k',linewidths=0.4)

        func_u = interp2d(x,y,u,kind='cubic')
        func_v = interp2d(x,y,v,kind='cubic')

        x = np.linspace(0,360,18)
        y = np.linspace(-90,90,16)

        u = func_u(x,y)
        v = func_v(x,y)

        u = np.flipud(u)
        v = np.flipud(v)

        # scale = 1/(0.001*np.amax(mag))

        ax.quiver(x,y,u,v,width=0.002,headwidth=3,scale=scale)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        c1 = plt.colorbar(p1, cax = cax, format='%.' +str(cformat) + 'f',orientation="horizontal")

        self.loadFieldColourOptions(c1,ax,cticks)
        self.loadFieldAxisOptions(ax)

    def PlotDissipation(self,extra_diss=[]):
        self.currentFig += 1

        norm = 100

        orbits = np.linspace(0,len(self.diss)+1,len(extra_diss))
        orbitst = np.linspace(0,len(self.diss),len(self.diss))

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        p1 = ax.plot(orbitst+0.5,self.diss,'k+',linewidth=0.8)
        p1 = ax.plot(orbitst+0.5,self.diss,'k--',linewidth=0.8,dashes=(5,5),label='Time Averaged Dissipation')

        if len(extra_diss)!=0:
            orbits_new = np.linspace(0,len(extra_diss)-1,len(extra_diss))

            plt.hold('on')
            p3 = ax.plot(orbits,extra_diss,color='r',linewidth=0.6,label='Instantaneous Dissipation')

        p1 = ax.plot(orbitst+0.5,self.diss,'k+',linewidth=0.8)
        p1 = ax.plot(orbitst+0.5,self.diss,'k--',linewidth=0.8,dashes=(5,5),label='Time Averaged Dissipation')

        # import scipy.fftpack
        # from scipy.interpolate import UnivariateSpline
        #
        # w = scipy.fftpack.rfft(extra_diss)
        # spectrum = w**2
        #
        # cutoff_idx = spectrum < (spectrum.max()/200)
        # w2 = w.copy()
        # w2[cutoff_idx] = 0
        # # w2[0] = 0
        # extra_diss = scipy.fftpack.irfft(w2)

        # sprange = np.linspace(0,len(self.diss),len(self.diss)*100)
        # sp = UnivariateSpline(orbits,extra_diss,s=4.5e-4)
        # sp = UnivariateSpline(orbitst,self.diss,s=8e-8)

        # plt.plot(sprange+0.5,sp(sprange))
        # plt.show()
        # plt.close()

        # p4 = ax.plot(sprange+0.5,sp(sprange),color='b',linewidth=0.6,label='Spline Fit')

        ax.set_xlim([min(orbits),max(orbits)])

        ax.set_xlabel('Orbit',fontsize=self.labelSize+1)
        ax.set_ylabel('Dissipation (\si{\watt\per\metre\squared})',fontsize=self.labelSize+1)

        ax.tick_params(axis='both', which='major', labelsize=self.slabelSize)
        ax.set_title('Enceladus Obliquity Tide Dissipation',fontsize=self.titleSize)

        legend = plt.legend()
        legend.get_frame().set_linewidth(0.5)
        for label in legend.get_texts():
            label.set_fontsize(10)
        for label in legend.get_lines():
            label.set_linewidth(0.8)


    def SetSuperTitle(self, title):
        self.fig.suptitle(title,fontsize=self.supTitleSize)

    def ShowFig(self):
        plt.show()
        plt.close()

    def SaveFig(self):
        self.fig.savefig(self.saveName +  '.pdf',bbox_inches='tight',pad_inches=1)
        self.fig.savefig(self.saveName +  '.png', format='PNG',bbox_inches='tight',dpi=300)

if __name__=="__main__":

    option = int(sys.argv[1])

    if option == 1:

        print("Option " + str(option) + "selected: Plotting displacement and velocity.\n")

        num = int(sys.argv[2])

        ODIS = ODISPlot(orbitnum=num,figdim=(1,2),figsize=(10,3),saveName="Test")

        ODIS.ReadFieldData("displacement")
        ODIS.PlotField(cticks = [])

        data_v = ODIS.ReadFieldData("north_vel")
        data_u = ODIS.ReadFieldData("east_vel")

        ODIS.PlotVelocity(data_u,data_v,cformat=2,scale=1/0.1,cticks=[])

        ODIS.SetSuperTitle("Test")

        ODIS.ShowFig()

        # ODIS.SaveFig()

    elif option == 2:

        print("Option " + str(option) + "selected: Plotting dissipation over time.\n")

        ODIS = ODISPlot(orbitnum=0,figdim=(1,1),figsize=(10,5),saveName="Dissipation")

        # n = ODIS.ReadDissipationData("/source/ODIS_dev/ODIS/Energy/kinetic_energy.txt")
        #
        # ODIS.ReadDissipationData("/source/ODIS_dev/ODIS/Energy/kinetic_energy_orb_avg.txt")

        n = ODIS.ReadDissipationData(os.getcwd() + "/Energy/diss_energy.txt")

        ODIS.ReadDissipationData(os.getcwd() + "/Energy/diss_energy_orb_avg.txt")

        #plt.hold('on')
        ODIS.PlotDissipation(extra_diss=n)

        ODIS.ShowFig()

        ODIS.SaveFig()

    elif option == 3:

        print("Option " + str(option) + "selected: Plotting velocity.\n")

        num = int(sys.argv[2])

        ODIS = ODISPlot(orbitnum=num,figdim=(1,2),figsize=(10,3),saveName="Test")

        data_v = ODIS.ReadFieldData("north_vel")
        ODIS.PlotField(cticks = [])
        data_u = ODIS.ReadFieldData("east_vel")
        ODIS.PlotField(cticks = [])

        ODIS.SetSuperTitle("Test")

        ODIS.ShowFig()

        # ODIS.SaveFig()

    elif option == 4:

        print("Option " + str(option) + "selected: Plotting displacement and dissipated energy.\n")

        num = int(sys.argv[2])

        ODIS = ODISPlot(orbitnum=num,figdim=(1,2),figsize=(10,3),saveName="Test")

        ODIS.ReadFieldData("displacement")
        ODIS.PlotField(cticks = [])

        data_d = ODIS.ReadFieldData("dissipation")

        ODIS.cmap = cm.magma

        ODIS.PlotField(data_d)

        ODIS.SetSuperTitle("Test")

        ODIS.ShowFig()


    else:

        print("Argument one must be either option 1-5.")

        #ODIS.SaveFig()

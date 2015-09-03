import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy as sc
from scipy.interpolate import interp2d
import os
import sys

class ODISPlot:
    def __init__(self,orbitnum=0,figdim=(1,2),figsize=(10,3),saveName="ODISPLOT"):
        # Set matplotlib text options
        plt.rc('font', family='serif')
        plt.rc('font',serif='Palatino')
        # for Palatino and other serif fonts use:
        # rc('font',**{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        plt.rc('text.latex',preamble='\\usepackage{siunitx}')

        # Set matplotlib rc options
        plt.rcParams['axes.linewidth'] = 0.8

        # Set colourmap
        self.cmap=cm.winter

        # Set label and title sizes
        self.labelSize = 10
        self.slabelSize = 9
        self.titleSize = 11
        self.supTitleSize = 12

        # Set plotting orbit
        self.orbit = orbitnum

        # Set figure dimensions
        self.currentFig = 0
        self.figSize = figsize
        self.figDim = figdim

        # Set figure save name
        self.saveName = saveName

        # Define figure
        self.fig = plt.figure(1,figsize=self.figSize)

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

        name = ["/EastVelocity/u_vel_"+orbitnum+".txt", "/NorthVelocity/v_vel_"+orbitnum+".txt", "/Displacement/eta_"+orbitnum+".txt"]

        if field=="displacement":
            name = name[2]
            self.caxisLabel = "Dispclacement, $\\eta$ (m)"
        elif field=="north_vel":
            name = name[1]
            self.caxisLabel = "North Velocity, $v$ (\\si{\\metre\\per\\second})"
        elif field=="east_vel":
            name = name[0]
            self.caxisLabel = "East Velocity, $v$ (\\si{\\metre\\per\\second})"
        else:
            print("Field " + field + "not recognised. Exiting.")
            sys.exit()

        path = directory + name
        self.data = []
        if os.path.exists(path):
            print("Retrieving data from:", path)
            init = open(path,'r')
            lines = init.readlines()
            init.close()

            for line in lines:
                line = line.split('\t')
                self.data.append([float(j) for j in line[:-1]])
            self.data = np.array(self.data)

        self.x = np.linspace(0,360,len(self.data[0]))
        self.y = np.linspace(90,-90,len(self.data))

        return self.data

    def ReadDissipationData(self,directory = os.getcwd()):
        self.diss = np.loadtxt(directory + "/diss.txt")[1:] # In Watts
        self.diss *= 1e3 # Convert to mWatts

    def PlotField(self,cformat=2,cticks=[]):
        self.currentFig += 1

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        self.loadFieldAxisOptions(ax)

        p1 = ax.pcolormesh(self.x,self.y,self.data,cmap=self.cmap,rasterized=False)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        c1 = plt.colorbar(p1, cax = cax, format='%.' +str(cformat) + 'f',orientation="horizontal")

        cs1 = ax.contour(self.x,self.y,self.data,8,colors='k',linewidths=0.4)

        self.loadFieldColourOptions(c1,ax,cticks)


    def PlotVelocity(self,u,v,cformat=3,scale=1/0.02,cticks=[]):
        self.currentFig += 1
        self.caxisLabel = "Velocity, $\\left|\\mathbf{u}\\right|$ (\\si{\\metre\\per\\second})"

        u=u[:-1][:]

        mag = np.sqrt(u**2+v**2)

        x = np.linspace(0,360,len(u[0]))
        y = np.linspace(-90,90,len(u))

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        p1 = ax.pcolormesh(x,y,mag,cmap=self.cmap,vmin=0,rasterized=False)
        cs1 = ax.contour(x,y,mag,8,colors='k',linewidths=0.4)

        func_u = interp2d(x,y,u,kind='cubic')
        func_v = interp2d(x,y,v,kind='cubic')

        x = np.linspace(0,360,18)
        y = np.linspace(-90,90,16)

        u = func_u(x,y)
        v = func_v(x,y)

        ax.quiver(x,y,u,v,pivot='middle',width=0.002,headwidth=3,scale=scale)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        c1 = plt.colorbar(p1, cax = cax, format='%.' +str(cformat) + 'f',orientation="horizontal")

        self.loadFieldColourOptions(c1,ax,cticks)
        self.loadFieldAxisOptions(ax)

    def PlotDissipation(self):
        self.currentFig += 1

        orbits = np.linspace(0,len(self.diss)-1,len(self.diss))
        print(orbits,self.diss)

        ax = self.fig.add_subplot(self.figDim[0],self.figDim[1],self.currentFig)
        p1 = ax.plot(orbits,self.diss,color='k',linewidth=0.6)
        p2 = ax.plot([orbits[0],orbits[-1]],[0.0547,0.0547],linestyle='--',color='k',linewidth=0.6)
        ax.set_yscale("log")
        ax.set_xlim([min(orbits),max(orbits)])

        ax.set_xlabel('Orbit',fontsize=self.labelSize+1)
        ax.set_ylabel('Dissipation (\si{\milli\watt\per\metre\squared})',fontsize=self.labelSize+1)

        ax.tick_params(axis='both', which='major', labelsize=self.slabelSize)
        ax.set_title('Obliquity Tide Dissipation',fontsize=self.titleSize)

        print(abs(0.0547-self.diss[-1])/0.0547*100)


    def SetSuperTitle(self, title):
        self.fig.suptitle(title,fontsize=self.supTitleSize)

    def ShowFig(self):
        plt.show()
        plt.close()

    def SaveFig(self):
        self.fig.savefig(self.saveName +  '.pdf',bbox_inches='tight')
        self.fig.savefig(self.saveName +  '.png', format='PNG',bbox_inches='tight',dpi=800)

if __name__=="__main__":

    option = int(sys.argv[1])

    if option == 1:

        print("Option 1 selected: Plotting displacement and velocity.\n")

        num = int(sys.argv[2])

        ODIS = ODISPlot(orbitnum=num,figdim=(1,2),figsize=(10,3),saveName="Test")

        ODIS.ReadFieldData("displacement")
        ODIS.PlotField(cticks = [])

        data_v = ODIS.ReadFieldData("north_vel")
        data_u = ODIS.ReadFieldData("east_vel")

        ODIS.PlotVelocity(data_u,data_v,cformat=2,scale=1/0.3,cticks=[])

        ODIS.SetSuperTitle("Test")

        ODIS.ShowFig()

    elif option == 2:

        ODIS = ODISPlot(orbitnum=0,figdim=(1,1),figsize=(4,3),saveName="obliq_diss")

        ODIS.ReadDissipationData()

        ODIS.PlotDissipation()

        ODIS.ShowFig()

        #ODIS.SaveFig()

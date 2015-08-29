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

    def ReadData(self,field="displacement",directory = os.getcwd()):
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

    def SetSuperTitle(self, title):
        self.fig.suptitle(title,fontsize=self.supTitleSize)

    def ShowFig(self):
        plt.show()
        plt.close()

    def SaveFig(self):
        self.fig.savefig(self.saveName +  '.pdf',bbox_inches='tight')
        self.fig.savefig(self.saveName +  '.png', format='PNG',bbox_inches='tight',dpi=800)

if __name__=="__main__":

    ODIS = ODISPlot(orbitnum=20,figdim=(1,2),figsize=(10,3))

    ODIS.ReadData("displacement")
    ODIS.PlotField(cticks = [-4,-2,0,2,4,6,8])

    data_v = ODIS.ReadData("north_vel")
    data_u = ODIS.ReadData("east_vel")

    ODIS.PlotVelocity(data_u,data_v,cformat=2,scale=1/0.3)

    ODIS.SetSuperTitle("c) Combined Eccentricity and Obliquity Tides")

    ODIS.ShowFig()
    ODIS.SaveFig()

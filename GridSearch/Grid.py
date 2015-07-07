import time
import os
import glob
import shutil
import numpy as np
import pprint
from Process import Process

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm

import scipy.interpolate
from scipy.interpolate import griddata

class Grid:
    def __init__(self, dim):
        self.dim = dim
        self.grid = [[0 for j in range(dim[0])] for i in range(dim[1])]
        self.max_proc = dim[0]*dim[1]

        self.running = []
        self.complete = []
        self.queue = []
        self.all = []

    def PopulateGrid(self, hrange, arange):
        self.hrange = hrange
        self.arange = arange
        self.omega = 4.559e-6
        self.height = np.logspace(hrange[0],hrange[1],self.dim[1])
        #self.alpha = omega/(2*np.logspace(arange[0],arange[1],self.dim[0]))
        #self.alpha = self.omega/(2*(10**np.linspace(arange[0],arange[1],self.dim[0])))
        self.alpha = np.logspace(arange[0],arange[1],self.dim[0])

        count = hcount = acount = 0
        for h in self.height:
            for a in self.alpha:
                self.grid[hcount][acount] = Process(count+1,h,a,(acount,hcount),self.dim)
                self.all.append(self.grid[hcount][acount])
                acount += 1
                count += 1

            acount = 0
            hcount += 1

        print(self.height)
        print(self.alpha)

    def RunProcess(self, p,parent=None, is_child=False):
        def IsComplete(process):
            if (os.path.exists(process.directory + "\\OUTPUT.txt")):
                logFile = open(process.directory + "\\OUTPUT.txt", 'r')
                string = logFile.readlines()[-4]
                logFile.close()

                if string == "Simulation complete.\n":
                    return True
                else:
                    return False
            else:
                return False

        def IsStart(process):
            return os.path.exists(process.directory + "\\OUTPUT.txt")

        destination = ["\\EastVelocity\\u_vel_", "\\NorthVelocity\\v_vel_", "\\Displacement\\eta_"]
        final = ["u_vel.txt", "v_vel.txt", "eta.txt"]

        if IsStart(p):
            if IsComplete(p):
                self.complete.append(p)
                print(p.id, "is already complete.")
                self.queue.remove(p)
                return 0
            else:

                if os.path.exists(p.directory + "\\InitialConditions"):
                    shutil.rmtree(p.directory + "\\InitialConditions")
                    os.makedirs(p.directory + "\\InitialConditions")

                for name in range(len(destination)):
                    all_files = []
                    for file in glob.glob(p.directory + destination[name] + "*.txt"):
                        all_files.append(file)

                    file_list = all_files.copy()
                    for i in range(len(file_list)):
                        file_list[i] = int((file_list[i].split('.')[-2]).split('_')[-1])

                    m = max(file_list)
                    max_index = [x for x, j in enumerate(file_list) if j == m]

                    shutil.copy(all_files[max_index[0]], p.directory + "\\InitialConditions")
                    print("Copying", all_files[max_index[0]], "to " + p.directory + "\\InitialConditions")
                    os.rename(p.directory + "\\InitialConditions\\" + all_files[max_index[0]].split('\\')[-1], p.directory + "\\InitialConditions\\" + final[name])

                shutil.rmtree(p.directory + "\\EastVelocity")
                os.makedirs(p.directory + "\\EastVelocity")
                shutil.rmtree(p.directory + "\\NorthVelocity")
                os.makedirs(p.directory + "\\NorthVelocity")
                shutil.rmtree(p.directory + "\\Displacement")
                os.makedirs(p.directory + "\\Displacement")


        else:

            if is_child and parent != None:
                par_dir = parent.directory

                os.makedirs(p.directory + "\\InitialConditions")

                for name in range(len(destination)):
                    all_files = []
                    for file in glob.glob(par_dir + destination[name] + "*.txt"):
                        all_files.append(file)

                    file_list = all_files.copy()
                    for i in range(len(file_list)):
                        file_list[i] = int((file_list[i].split('.')[-2]).split('_')[-1])

                    m = max(file_list)
                    max_index = [x for x, j in enumerate(file_list) if j == m]

                    shutil.copy(all_files[max_index[0]], p.directory + "\\InitialConditions")
                    print("Copying", all_files[max_index[0]], "to " + p.directory + "\\InitialConditions")
                    os.rename(p.directory + "\\InitialConditions\\" + all_files[max_index[0]].split('\\')[-1], p.directory + "\\InitialConditions\\" + final[name])

        p.Run()
        self.running.append(p)
        print(p.id, "is now running.")
        self.queue.remove(p)
        return 1

    def CheckRunning(self):
        for proc in self.running:
            if proc.IsRun():
                self.running.remove(proc)
                print(proc.id, " finished.")
                self.complete.append(proc)

    def SolveGrid(self,total):
        hlen = len(self.grid)
        alen = len(self.grid[0])

        diagNum = alen + hlen - 1
        diagList = []

        x = 0
        ypos = 1
        for d in range(diagNum):
            y = 0
            x = d
            if d < alen:
                while (x >= 0) and (y < hlen):
                    diagList.append((x,y))
                    x -= 1
                    y += 1
            else:
                x = alen -1
                y = ypos
                while (y < hlen):
                    diagList.append((x,y))
                    x -= 1
                    y += 1
                ypos +=1

        #Create ID queue
        for pos in diagList:
            for h in range(hlen):
                for a in range(alen):
                    if pos == self.grid[h][a].node:
                        self.queue.append(self.grid[h][a])

        #Solve first node - others spawn from this
        if self.RunProcess(self.grid[0][0]) != 0: #not already run
            #Wait for process 1 to finish before entering main loop
            exit_code = self.running[0].Wait()

        while len(self.running)!= 0 or len(self.queue) != 0:
            #check complete list

            looped = 0

            while looped <= 8:
                if len(self.complete) > 0:
                    for p in self.complete:
                        if p.max_children != 0:

                            found_children = []

                            for child in self.queue:
                                if len(found_children)<p.max_children:
                                    if p.IsChild(child,alen):
                                        #child process found
                                        found_children.append(child)
                                else:
                                    break


                            if len(self.running) < total and len(found_children) > 0:
                                for newp in found_children:
                                    if not newp.IsRun() and len(self.running) < total:
                                        self.RunProcess(newp,parent=p,is_child=True, )

                    self.CheckRunning()

                else:
                    self.CheckRunning()

                looped+=1
                time.sleep(15)

            print("Queued processes: ", self.queue)
            print("Running processes: ", self.running)
            print("Completed processes: ", self.complete, "\n")

    def CollectResults(self):
        all_dirs = []

        results = open(os.getcwd()+"\\grid_results.txt",'w')

        self.diss = []
        for p in self.all:
            if os.path.isfile(p.results_dir + "\\diss.txt"):
                diss_file = open(p.results_dir + "\\diss.txt")
                lines = diss_file.readlines()
                results.write(str(p.h) + "\t" + str(p.a) + "\t" + lines[-1])
                self.diss.append(float(lines[-1]))
                diss_file.close()

        results.close()

        self.PlotResults()

    def PlotResults(self):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import matplotlib

        h,a,diss = np.loadtxt(os.getcwd()+"\\grid_results.txt").T #Transposed for easier unpacking
        diss = diss*1e-3

        h_dat = np.zeros(len(self.height))
        a_dat = np.zeros(len(self.alpha))

        h_dat[0] = h[0]
        count = 0
        for i in range(len(h)):
            if h_dat[count] != h[i]:
                count += 1
                h_dat[count] = h[i]

        a_dat[0] = a[0]
        count = 0
        for i in range(len(a)):
            add = True
            for j in range(len(a_dat)):
                if a_dat[j] == a[i]:
                    add = False

            if add:
                count += 1
                a_dat[count] = a[i]

        Z = np.zeros((len(a_dat),len(h_dat)))
        count = 0
        for j in range(len(h_dat)):
            for i in range(len(a_dat)):
                Z[i][j] = diss[count]
                count += 1


        res=400
        y = np.logspace(self.hrange[0],self.hrange[1],res)
        #self.alpha = omega/(2*np.logspace(arange[0],arange[1],self.dim[0]))
        x = self.omega/(2*(10**np.linspace(self.arange[0],self.arange[1],res)))

        X2, Y2 = np.meshgrid(a_dat,h_dat)
        X, Y = np.meshgrid(x,y)
        #diss2 = np.log10(diss)
        diss2 = np.log10(diss)
        #zi = griddata((a, h), diss2, (X, Y),method='linear')

        zi = griddata((a, h), diss2, (X, Y),method='cubic')

        fig = plt.figure(1,figsize=(7, 10), dpi=600)


        plt.rc('font', family='serif')
        plt.rc('font',serif='Palatino')
        # for Palatino and other serif fonts use:
        # rc('font',**{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)


        titleSize = 14
        titleSize2 = 12
        labelSize = 10
        tickSize = 10

        fig.suptitle('Ocean Dissipation for Eccentricity-Libration Tide', fontsize=titleSize)

        ax1 = fig.add_subplot(2,1,1)
        m1 = ax1.pcolormesh(X,Y,zi)
        c1 = plt.colorbar(m1, ax = ax1, format='%.1f')
        c1.set_label("$\\log_{10}$ (Dissipated Energy), (mW m$\displaystyle{^{-2}}$)",fontsize=labelSize)
        c1.outline.set_linewidth(0.7)
        c1.set_ticks([-8.0, -7,-6,-5, -4,-3])
        ax1.scatter(a,h,marker="+",s=8,color='k',alpha=0.7,linewidths=0.6)


        ax2 = fig.add_subplot(2,1,2)
        m2 = ax2.pcolor(X2,Y2,np.log10(Z.T))
        c2 = plt.colorbar(m2, ax = ax2, format='%.1f')
        c2.set_label("$\\log_{10}$ (Dissipated Energy), (mW m$\displaystyle{^{-2}}$)",fontsize=labelSize)
        c2.outline.set_linewidth(0.7)
        c2.set_ticks([-8.0, -7,-6,-5, -4])
        ax2.scatter(a,h,marker="+",s=8,color='k',alpha=0.7,linewidths=0.6)

        hmin = 10
        hmax = 500

        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax2.set_xscale("log")
        ax1.invert_xaxis()
        ax1.axis([a_dat.max(), a_dat.min(), hmin, hmax])
        ax2.axis([a_dat.max(), a_dat.min(), hmin, hmax])
        ax2.set_xlabel('$\\alpha$ (s$\displaystyle{^{-1}}$)',fontsize=labelSize)
        ax1.set_ylabel('Ocean Depth (m)',fontsize=labelSize)
        ax2.set_ylabel('Ocean Depth (m)',fontsize=labelSize)
        ax1.set_title('Interpolated Data',fontsize=titleSize2)
        ax2.set_title('Discrete Data',fontsize=titleSize2)

        yticks = [10,100,500]
        ax1.yaxis.set_ticks( yticks )
        ax1.yaxis.set_ticklabels( ['%d' % i for i in yticks] )


        ax2.yaxis.set_ticks( yticks )
        ax2.yaxis.set_ticklabels( ['%d' % i for i in yticks] )

        ax1.tick_params(axis='both', which='major', labelsize=tickSize)
        ax2.tick_params(axis='both', which='major', labelsize=tickSize)
        c1.ax.tick_params(axis='both', which='major', labelsize=tickSize)
        c2.ax.tick_params(axis='both', which='major', labelsize=tickSize)

        #plt.show()



        fig.savefig('Eccentricity-Libration.png', format='PNG',dpi=1200)

        #diss = diss*1e-3

        plt.close()

        fig2 = plt.figure(2, dpi=600)
        ax12 = fig2.add_subplot(1,1,1)

        labels = ['%.2e' % a_dat[i] for i in range(len(a_dat))]
        diss_interp = 10**zi[:,0]
        p1 = ax12.plot(h_dat,diss[:-1:8]*1e3,marker='+',label='$\\alpha = $ ' + labels[0] + ' s$^{-1}$')
        p2 = ax12.plot(h_dat,diss[3:-1:8]*1e3,marker='+',label='$\\alpha = $ ' + labels[3] + ' s$^{-1}$')
        p3 = ax12.plot(h_dat,diss[6:-1:8]*1e3,marker='+',label='$\\alpha = $ ' + labels[6] + ' s$^{-1}$')
        p4 = ax12.plot(y,(10**zi[:,0])*1e3,ls=':',color='k',label='$\\alpha = $ ' + labels[0] + ' s$^{-1}$: interpolated')
        ax12.set_xscale("log")
        ax12.set_yscale("log")
        ax12.set_ylabel("Dissipated Energy, (W m$^{-2}$)",fontsize=labelSize)
        ax12.set_xlabel("Ocean Depth, (m)",fontsize=labelSize)
        ax12.set_title('Ocean Dissipation for Eccentricity-Libration Tide',fontsize=titleSize)
        plt.ylim([1e-6, 1e1])
        plt.xlim([1e1, 600])
        plt.legend()
        #p2 = ax12.plot(y,zi[:,0]*1e3)
        #plt.show()
        fig2.savefig('Eccentricity-Libration_depth.png', format='PNG',dpi=1200)
        #pp = PdfPages('Eccentricity-Libration.pdf')
        #pp.savefig(fig)
        #pp.close()

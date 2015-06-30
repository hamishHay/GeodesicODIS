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
        omega = 4.559e-6
        self.height = np.logspace(hrange[0],hrange[1],self.dim[1])
        #self.alpha = omega/(2*np.logspace(arange[0],arange[1],self.dim[0]))
        self.alpha = omega/(2*(10**np.linspace(arange[0],arange[1],self.dim[0])))

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
        if is_child and parent != None:
            par_dir = parent.directory
            os.makedirs(p.directory + "\\InitialConditions")

            destination = ["\\EastVelocity\\u_vel_", "\\NorthVelocity\\v_vel_", "\\Displacement\\eta_"]
            final = ["u_vel.txt", "v_vel.txt", "eta.txt"]

            for name in range(len(destination)):
                all_files = []
                for file in glob.glob(par_dir + destination[name] + "*.txt"):
                    all_files.append(file)

                file_list = all_files.copy()
                for i in range(len(file_list)):
                    file_list[i] = (file_list[i].split('.')[-2]).split('_')[-1]

                m = max(file_list)
                max_index = [x for x, j in enumerate(file_list) if j == m]

                shutil.copy(all_files[max_index[0]], p.directory + "\\InitialConditions")
                os.rename(p.directory + "\\InitialConditions\\" + all_files[max_index[0]].split('\\')[-1], p.directory + "\\InitialConditions\\" + final[name])

        p.Run()
        self.running.append(p)
        print(p.id, "is now running.")
        self.queue.remove(p)

    def CheckRunning(self):
        for proc in self.running:
            if proc.IsRun():
                self.running.remove(proc)
                print(proc.id, " finished.")
                self.complete.append(proc)

    def SolveGrid(self,total):
        hlen = len(self.grid[0])
        alen = len(self.grid)

        diagNum = alen + hlen - 1
        diagList = []
        for d in range(diagNum):
            if d < alen:
                for i in range(d+1):
                    if i < hlen:
                        diagList.append((d-i,i))
            else:
                for i in range(d-hlen, hlen):
                    diagList.append((d-i,i))

        #Create ID queue
        for pos in diagList:
            for h in range(hlen):
                for a in range(alen):
                    if pos == self.grid[a][h].node:
                        self.queue.append(self.grid[a][h])

        #Solve first node - others spawn from this
        self.RunProcess(self.grid[0][0])

        #Wait for process 1 to finish before entering main loop
        exit_code = self.running[0].Wait()

        while len(self.complete)!=self.max_proc:
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

        #self.PlotResults()

    def PlotResults(self):

        h,a,diss = np.loadtxt(os.getcwd()+"\\grid_results.txt").T #Transposed for easier unpacking
        diss = diss*1e3
        print(diss)
        #nrows, ncols = #len(self.height), len(self.alpha)
        #data = diss.reshape((8, 8))

        #x, y = np.meshgrid(self.alpha, self.height)

        xi = a#np.linspace(self.alpha[0],self.alpha[-1],1000)
        yi = h#np.linspace(self.height[0],self.height[-1],1000)

        x = a#self.alpha.copy()
        y = h#self.height.copy()

        zi = griddata(x,y,diss,xi,yi,interp='linear')

        #fig = plt.figure()
        #ax = fig.add_subplot(1,1,1)


        plt.pcolor(xi, yi, zi, norm=LogNorm(vmin=0.003, vmax=zi.max()))
        #CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
        #                  vmax=abs(zi).max(), vmin=-abs(zi).max())
        plt.yscale("log")
        plt.xscale("log")
        #ax.set_yscale('log')
        plt.colorbar()  # draw colorbar
        # plot data points.
        plt.scatter(x, y, marker='o', c='b', s=5, zorder=10)
        #plt.xlim(x[0], x[-1])
        #plt.ylim(y[0], y[-1])
        #plt.title('griddata test (%d points)' % npts)
        plt.show()

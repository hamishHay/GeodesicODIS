import time
import os
import sys
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

    def IsComplete(self, process):
        if (os.path.exists(process.directory + "/OUTPUT.txt")):
            logFile = open(process.directory + "/OUTPUT.txt", 'r')
            string = logFile.readlines()[-4]
            logFile.close()
            if string == "Simulation complete.\n":
                return True
            else:
                return False
        else:
            return False

    def IsStart(self, process):
        return os.path.exists(process.directory + "/Displacement/eta_0.txt")

    def MakeDir(self,directory):
        if os.path.exists(directory):
            print(directory,"already exists. Proceeding to overwrite.")
            shutil.rmtree(directory)

        os.makedirs(directory)

    def CopyInitialConditions(self,copy_dir,place_dir):
        field_dir = ["/EastVelocity/u_vel_", "/NorthVelocity/v_vel_", "/Displacement/eta_"]
        field_txt = ["u_vel.txt", "v_vel.txt", "eta.txt"]

        ### For each field, find latest output file and copy into InitialConditions
        for k in range(len(field_dir)):

            ### List all files in field_dir[k]
            all_files = []
            for file in glob.glob(copy_dir + field_dir[k] + "*.txt"):
                all_files.append(file)

            ### Strip file names to get output number and store in file_list
            file_list = all_files.copy()
            for i in range(len(file_list)):
                file_list[i] = int((file_list[i].split('.')[-2]).split('_')[-1])

            ### If no files in folder, then process has not started
            if len(file_list) == 0:
                print("No initial conditions found for", place_dir)
                return 0

            ### Files found, find the index of the highest output number
            else:
                m = max(file_list)
                max_index = [x for x, j in enumerate(file_list) if j == m]
                print("Initial conditions found for", place_dir)

            ### Copy initial condition from field_dir into initial conditions, and rename
            shutil.copy(all_files[max_index[0]], place_dir + "/InitialConditions")
            print("Copying", all_files[max_index[0]], "to " + place_dir + "/InitialConditions")
            os.rename(place_dir + "/InitialConditions/" + all_files[max_index[0]].split('/')[-1], place_dir + "/InitialConditions/" + field_txt[k])

    def RunProcess(self, p,parent=None, is_child=False):
        ### To run a process, three steps are taken:
        ### 1. Check if the process is already started.
        ###    This is true if the process directory contains an eta_0.txt file.
        ### 2. Check if the process is finished.
        ###    Checks for the simulation complete tag at the end of OUTPUT.txt (no very robust)
        ### 3. Run process if not complete.
        ###    This involes checking for initial conditions to use as the simulation starting point
        ###    If it is started, then copy last output from field_dir into initial condition folder.
        ###    If it is not started, then do the above but using the par_dir instead.

        ### Check for signs of started process
        if self.IsStart(p):
            ### If started, then check for signs of completeness
            if self.IsComplete(p) and self.refine == False:
                self.complete.append(p)
                print(p.id, "is already complete: "+ str(p.h) +", " + str(p.a))
                self.queue.remove(p)
                return 0

            ### Started, but not complete - copy initial conditions from own directory
            else:
                ### Remove current inital conditions as they are now obsolete
                self.MakeDir(p.directory + "/InitialConditions")

                ### Copy initial conditions in process to process
                self.CopyInitialConditions(p.directory,p.directory)

                ### Process is now ready to run, so remove old values
                shutil.rmtree(p.directory + "/EastVelocity")
                shutil.rmtree(p.directory + "/NorthVelocity")
                shutil.rmtree(p.directory + "/Displacement")

        ### Process is not started - find initial conditions from parent
        else:
            ### Check to see if the process is a child and that the child has a parent
            if is_child and parent != None:

                ### Create initial condition directory
                self.MakeDir(p.directory + "/InitialConditions")

                ### Copy initial conditions in parent to process
                self.CopyInitialConditions(parent.directory, p.directory)

        if p.id != 1:
            self.RemapInitCondition(p)

        p.Run()
        self.running.append(p)
        print(p.id, "is now running. h:", p.h, "a:", p.a)
        self.queue.remove(p)
        return 1

    def CheckRunning(self):
        for proc in self.running:
            if proc.IsRun():
                self.running.remove(proc)
                print(proc.id, " finished.")
                self.complete.append(proc)

    def SolveGrid(self,total,refine=False):
        hlen = len(self.grid)
        alen = len(self.grid[0])
        self.refine = refine

        #Create ID queue
        for h in range(hlen):
            for a in range(alen):
                self.queue.append(self.grid[h][a])

        #Assess Completed
        for process in self.queue:
            if self.IsComplete(process):
                self.queue.remove(process)
                self.complete.append(process)

        #Solve first node - others spawn from this
        if not self.IsComplete(self.grid[0][0]):
            if self.RunProcess(self.grid[0][0]) != 0 and self.refine == False: #not already run
                #Wait for process 1 to finish before entering main loop
                exit_code = self.running[0].Wait()

        self.rerun = False
        while len(self.running)!= 0 or len(self.queue) != 0:
            #check complete list

            looped = 0

            while looped <= 120:
                # If refining Grid, then no need to run when parent complete
                if self.refine or self.rerun:
                    for p in self.queue:
                        if len(self.running) < total and len(self.queue) > 0:
                            self.RunProcess(p)

                    self.CheckRunning()

                else:
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

                                            #Attempt to run process
                                            self.RunProcess(newp,parent=p,is_child=True)

                        self.CheckRunning()

                    else:
                        self.CheckRunning()

                looped+=1
                time.sleep(2)

            print("Queued processes: ", self.queue)
            print("Running processes: ", self.running)
            print("Completed processes: ", self.complete)
            print("Number remaining: ", len(self.queue) + len(self.running), "\n")

    def ResetRange(self, a1, a2, h1, h2):
        self.rerun = True
        for p in self.all:
            if (p.a > a1) and (p.a < a2) and (p.h > h1) and (p.h < h2):
                rfile = open(p.directory+"/OUTPUT.txt",'r')
                lines = rfile.readlines()
                lines = lines[:-5] #remove last five lines
                rfile.close()

                wfile = open(p.directory+"/OUTPUT.txt",'w')
                wfile.writelines([item for item in lines])
                wfile.close()

    def DeleteRange(self, a1, a2, h1, h2):
        import shutil
        self.rerun = True
        for p in self.all:
            if (p.a > a1) and (p.a < a2) and (p.h > h1) and (p.h < h2):
                try:
                    shutil.rmtree(p.directory)
                except WindowsError:
                    print("he")

    def MakeComplete(self, a1, a2, h1, h2):
        #self.rerun = True
        for p in self.all:
            if (p.a > a1) and (p.a < a2) and (p.h > h1) and (p.h < h2):
                rfile = open(p.directory+"/OUTPUT.txt",'r')
                lines = rfile.readlines()
                lines.append('\n')
                lines.append('Simulation complete.\n')
                lines.append('\n')
                lines.append('End time: 		made complete\n')
                lines.append('Total iterations: 	made complete\n')
                rfile.close()

                wfile = open(p.directory+"/OUTPUT.txt",'w')
                wfile.writelines([item for item in lines])
                wfile.close()


    def CollectResults(self):
        from scipy import stats as st

        def findMid(arr):
            first=0
            second=0

            slopeSign = (arr[-1]-arr[-2])/abs((arr[-1]-arr[-2]))

            for i in range(len(arr)-1,1,-1):
                nSlopeSign = (arr[i]-arr[i-1])/abs((arr[i]-arr[i-1]))
                if nSlopeSign != slopeSign:
                    first = i
                    slopeSign = nSlopeSign
                    break

            for i in range(first,1,-1):
                nSlopeSign = (arr[i]-arr[i-1])/abs((arr[i]-arr[i-1]))
                if nSlopeSign != slopeSign:
                    second = i
                    break

            return int((first-second)/2)

        all_dirs = []

        results = open(os.getcwd()+"/grid_results.txt",'w')

        self.diss = []
        for p in self.all:
            if os.path.isfile(p.results_dir + "/diss.txt"):
                diss_file = open(p.results_dir + "/diss.txt")
                lines = diss_file.readlines()

                #print(p.results_dir)
                flines = [float(x) for x in lines]
                flines = np.array(flines[1:])
                bot = 30

                resid = abs(flines[1:-1] - flines[0:-2])

                converged = True
                if (len(flines)>3):
                    for i in range(len(resid)-1,len(resid)-3,-1):
                        #print(resid[i],i)
                        if resid[i] > 5e-9:
                            converged=False
                            break
                else:
                    converged=False

                if not converged:# and len(flines)>=12:#p.a <1e-7:# and np.log10(flines[-1])<-4.3:
                    bot=1
                    x = np.linspace(bot,len(flines),len(flines[bot:]))

                    m, c, r, pval, e = st.linregress(x,flines[bot:])
                    #print(p.h, p.a, len(flines), m)

                    if abs(m)>=1e-6:
                        val = flines[-1]
                    elif abs(m) < 1e-6 and len(flines)>20:
                        #mid = findMid(flines)
                        mid = len(flines)/2
                        #val = flines[mid]
                        #val = flines[-1]
                        val = np.median(flines[mid:])
                    elif len(flines)>2:
                        val = flines[-1]
                        #val = np.median(flines[bot:])
                    #else:
                    #    val = flines[-1]


                else:
                    val = flines[-1]
                #print(val)
                results.write(str(p.h) + "\t" + str(p.a) + "\t" + str(val) + "\n")#lines[-1])
                self.diss.append(str(val))
                diss_file.close()

        results.close()

        #self.PlotResults()

    def PlotResults(self):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import matplotlib
        from matplotlib import cm

        def GetVal(name="/grid_results.txt"):
            self.h,self.a,diss = np.loadtxt(os.getcwd()+name).T #Transposed for easier unpacking
            self.diss = diss*1e-3
            print(self.diss)
            #data = np.zeros((len(a_dat),len(h_dat)))

            #h_dat = np.zeros(len(self.height))
            #a_dat = np.zeros(len(self.alpha))

            #data = np.zeros((len(a_dat),len(h_dat)))

            h_dat[0] = self.h[0]
            count = 0
            for i in range(len(self.h)):
                if h_dat[count] != self.h[i]:
                    count += 1
                    h_dat[count] = self.h[i]

            a_dat[0] = self.a[0]
            count = 0
            for i in range(len(self.a)):
                add = True
                for j in range(len(a_dat)):
                    if a_dat[j] == self.a[i]:
                        add = False

                if add:
                    count += 1
                    a_dat[count] = self.a[i]

            data = np.zeros((len(a_dat),len(h_dat)))
            count = 0
            for j in range(len(h_dat)):
                for i in range(len(a_dat)):
                    #data[i][j] = self.diss[count]*1e3
                    # try:
                    #     data[i][j] = self.diss[count]*1e3
                    # except IndexError:
                    #     data[i][j] = 1
                    #     #print("Is NaN: self.h: " + h_dat[j] + ", a: " + a_dat[i])

                    count += 1

            return data


        h_dat = np.zeros(len(self.height))
        a_dat = np.zeros(len(self.alpha))

       # Z2 = GetVal(name="/grid_results1.txt")

        Z = GetVal(name="/grid_results.txt")
        print(Z)
        # while True in ps.isnull(Z):# or 0.0 in Z:
        #     print(ps.isnull(data))
        #     time.sleep(2)
        #     Z = GetVal(name="/grid_results.txt")
        h = self.h
        a = self.a

        #Z = Z-Z2
        cmap = cmap=cm.gnuplot2
        #print(self.diss)
        fig2 = plt.figure()
        ax12 = fig2.add_subplot(1,1,1)
        pl = ax12.scatter(a,h,c=np.log10(self.diss*1e3),cmap=cmap,lw = 0.5)

        #pl = ax12.plot(a, np.log10(self.diss*1e3),marker='+')

        c1 = plt.colorbar(pl, ax = ax12, format='%.1f')
        #ax12.set_clim(vmin=-5.5, vmax=0.4)
        ax12.set_xscale("log")
        ax12.set_yscale("log")
        plt.ylim([10**self.hrange[0],10**self.hrange[1]])
        plt.xlim([10**self.arange[1],10**self.arange[0]])
        plt.show()

        sys.exit()
        res=400
        y = np.logspace(self.hrange[0],self.hrange[1],res)
        #self.alpha = omega/(2*np.logspace(arange[0],arange[1],self.dim[0]))
        x = np.logspace(self.arange[0],self.arange[1],res)

        X2, Y2 = np.meshgrid(a_dat,h_dat)
        X, Y = np.meshgrid(x,y)
        #diss2 = np.log10(diss)
        diss2 = np.log10(self.diss)
        #zi = griddata((a, h), diss2, (X, Y),method='linear')

        #zi = griddata((a, h), diss2, (X, Y),method='cubic')

        import scipy as sc
#np.fliplr(np.log10(Z.T))
        #Much better than griddata
        func = sc.interpolate.interp2d(a_dat,h_dat,np.fliplr(np.log10(Z.T)),kind='linear')
        #func = sc.interpolate.SmoothBivariateSpline(np.sort(a),np.sort(h),diss2)
        zi = func(x,y)

        plt.rc('font', family='serif')
        plt.rc('font',serif='Palatino')
        # for Palatino and other serif fonts use:
        # rc('font',**{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)

        name = 'Full'


        titleSize = 14
        titleSize2 = 12
        labelSize = 10
        tickSize = 10

        fig = plt.figure(1,figsize=(7, 10), dpi=80)
        fig.suptitle('Ocean Dissipation for ' + name + ' Tide', fontsize=titleSize)
        ticks = [-10.0,-9.0, -8.0, -7,-6,-5, -4,-3,-2,-2,-1,0,1]

        ax1 = fig.add_subplot(2,1,1)
        m1 = ax1.pcolormesh(X,Y,zi,cmap=cmap)
        #ax1.contour(X,Y,zi,10,cmap=cm.gray,linewidths=0.8)
        c1 = plt.colorbar(m1, ax = ax1, format='%.1f')
        c1.set_label("$/log_{10}$ (Dissipated Energy), (W m$\displaystyle{^{-2}}$)",fontsize=labelSize)
        c1.outline.set_linewidth(0.7)
        c1.set_ticks(ticks)
        ax1.scatter(a,h,marker="+",s=7,color='k',alpha=0.7,linewidths=0.4)

        # for j in range(len(h_dat)):
        #     for i in range(len(a_dat)):
        #         if Z[i][j] != Z2[i][j]:
        #             ax1.scatter(a_dat[i],h_dat[j],marker="+",s=8,color='w',alpha=0.7,linewidths=0.6)
        #         else:
        #             ax1.scatter(a_dat[i],h_dat[j],marker="+",s=8,color='k',alpha=0.7,linewidths=0.6)


        ax2 = fig.add_subplot(2,1,2)
        m2 = ax2.pcolor(X2,Y2,np.log10(Z.T),cmap=cmap)
        c2 = plt.colorbar(m2, ax = ax2, format='%.1f')
        c2.set_label("$/log_{10}$ (Dissipated Energy), (W m$\displaystyle{^{-2}}$)",fontsize=labelSize)
        c2.outline.set_linewidth(0.7)
        c2.set_ticks(ticks)
        ax2.scatter(a,h,marker="+",s=7,color='k',alpha=0.7,linewidths=0.4)

        # for j in range(len(h_dat)):
        #     for i in range(len(a_dat)):
        #         if Z[i][j] != Z2[i][j]:
        #             ax2.scatter(a_dat[i],h_dat[j],marker="+",s=8,color='w',alpha=0.7,linewidths=0.6)
        #         else:
        #             ax2.scatter(a_dat[i],h_dat[j],marker="+",s=8,color='k',alpha=0.7,linewidths=0.6)

        hmin = 1
        hmax = 10000

        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax2.set_xscale("log")
        #ax1.invert_xaxis()
        ax1.axis([a_dat.min(), a_dat.max(), hmin, hmax])
        ax2.axis([a_dat.min(), a_dat.max(), hmin, hmax])
        ax2.set_xlabel('$/alpha$ (s$\displaystyle{^{-1}}$)',fontsize=labelSize)
        ax1.set_ylabel('Ocean Depth (m)',fontsize=labelSize)
        ax2.set_ylabel('Ocean Depth (m)',fontsize=labelSize)
        ax1.set_title('Interpolated Data',fontsize=titleSize2)
        ax2.set_title('Discrete Data',fontsize=titleSize2)

        yticks = [1,10,100,1000,10000]
        ax1.yaxis.set_ticks( yticks )
        ax1.yaxis.set_ticklabels( ['%d' % i for i in yticks] )


        ax2.yaxis.set_ticks( yticks )
        ax2.yaxis.set_ticklabels( ['%d' % i for i in yticks] )

        ax1.tick_params(axis='both', which='major', labelsize=tickSize)
        ax2.tick_params(axis='both', which='major', labelsize=tickSize)
        c1.ax.tick_params(axis='both', which='major', labelsize=tickSize)
        c2.ax.tick_params(axis='both', which='major', labelsize=tickSize)

        plt.show()

        # fig.savefig(name + '.pdf', format='PDF')
        fig.savefig(name + '.png', format='PNG',dpi=1200)

        #self.diss = diss*1e-3

        plt.close()
#
#         fig2 = plt.figure(2, dpi=80)
#         ax12 = fig2.add_subplot(1,1,1)
#
#         labels = ['%.2e' % a_dat[i] for i in range(len(a_dat))]
#         diss_interp = 10**zi[:,0]
#         p1 = ax12.plot(h_dat,Z[6],marker='+',label='$/alpha = $ ' + labels[6] + ' s$^{-1}$')
#         p2 = ax12.plot(h_dat,Z[7],marker='+',label='$/alpha = $ ' + labels[7] + ' s$^{-1}$')
#         p1 = ax12.plot(h_dat,Z[8],marker='+',label='$/alpha = $ ' + labels[8] + ' s$^{-1}$')
#         p2 = ax12.plot(h_dat,Z[9],marker='+',label='$\\alpha = $ ' + labels[9] + ' s$^{-1}$')
#         p3 = ax12.plot(h_dat,Z[10],marker='+',label='$\\alpha = $ ' + labels[10] + ' s$^{-1}$')
#         p3 = ax12.plot(h_dat,Z[11],marker='+',label='$\\alpha = $ ' + labels[11] + ' s$^{-1}$')
#         p3 = ax12.plot(h_dat,Z[12],marker='+',label='$\\alpha = $ ' + labels[12] + ' s$^{-1}$')
#         #p4 = ax12.plot(h_dat,Z[9],marker='+',label='$\\alpha = $ ' + labels[9] + ' s$^{-1}$')
#         #p4 = ax12.plot(y,(10**zi[:,5]),ls=':',color='k',label='$\\alpha = $ ' + labels[5] + ' s$^{-1}$: interpolated')
#         ax12.set_xscale("log")
#         ax12.set_yscale("log")
#         ax12.set_ylabel("Dissipated Energy, (W m$^{-2}$)",fontsize=labelSize)
#         ax12.set_xlabel("Ocean Depth, (m)",fontsize=labelSize)
#         ax12.set_title('Ocean Dissipation for ' + name + ' Tide',fontsize=titleSize)
#         plt.ylim([1e-8, 1e1])
#         plt.xlim([1, 10000])
#         plt.legend()
#         #p2 = ax12.plot(y,zi[:,0]*1e3)
#         plt.show()
#         # fig2.savefig(name + '_depth.pdf', format='PDF')
#         # fig2.savefig(name + '_depth.png', format='PNG',dpi=1200)
#
#         #pp = PdfPages('Eccentricity-Libration.pdf')
#         #pp.savefig(fig)
#         #pp.close()

    def RemapInitCondition(self,p):
        import scipy as sc

        name = ["/u_vel.txt", "/v_vel.txt", "/eta.txt"]
        for i in range(3):
            path = p.directory + "/InitialConditions/" + name[i]
            #path = os.getcwd() + name[i]
            data = []
            #p.directory + "/InitialConditions"
            if os.path.exists(path):
                init = open(path,'r')
                lines = init.readlines()
                init.close()

                for line in lines:
                    line = line.split('\t')
                    data.append([float(j) for j in line[:-1]])
                data = np.array(data)

            lat = np.linspace(0,len(data),len(data))
            lon = np.linspace(0,len(data[0]),len(data[0]))


            if i == 1:
                latnew = np.linspace(0,len(data),p.latnum)
            else:
                latnew = np.linspace(0,len(data),p.latnum+1)
            lonnew = np.linspace(0,len(data[0]),p.lonnum)

            func = sc.interpolate.interp2d(lon,lat,data,kind='cubic')
            #func = sc.interpolate.SmoothBivariateSpline(np.sort(a),np.sort(h),diss2)
            zi = func(lonnew,latnew)

            init = open(path,'w')
            for line in zi:
                for val in line:
                    init.write(str(val) + '\t')
                init.write('\n')

            init.close()

#! python

import subprocess as sub
import time
import os
import shutil
import random as rand
import numpy as np
import pprint

exect = "C:\\Users\\Hamish\\Documents\\GitHub\\GridTest\\runProcess.py"

class Process:
    def __init__(self, ID, height, alpha, Node, dim):
        self.id = ID
        self.h = height
        self.a = alpha
        self.start = False
        self.complete = False
        self.node = Node
        self.childrenPos = []

        if self.node[0] == 0 and self.node[1] == dim[1] - 1: #top and left
            self.childrenPos.append((self.node[0]+1,self.node[1]))
            self.max_children = 1
        elif self.node[0] == 0: #left
            self.childrenPos.append((self.node[0]+1,self.node[1]))
            self.childrenPos.append((self.node[0],self.node[1]+1))
            self.max_children = 2
        elif self.node[0] == dim[0] - 1: #right
            self.max_children = 0
        else: #interior
            self.childrenPos.append((self.node[0]+1,self.node[1]))
            self.max_children = 1
            
    def __repr__(self):
        return str(self.id)#"ID="+str(self.id)+",h="+str(self.h)+",a=" + str(self.a)+",comp="+str(self.complete)

    def __str__(self):
        return str(self.id)

    def CreateDir(self):
        self.directory = os.getcwd() + "\\h" + str(self.h) + "_alpha" + str(self.a)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            shutil.copy("ODIS.exe", self.directory)
            #shutil.copy("input.in", self.directory)
            self.CreateInputFile()

    def CreateInputFile(self):
        f = open(self.directory+"\\input.in",'w')
        f.write("ocean thickness; \t \t \t " + str(self.h) +";	\t \t \t h; \n")
        f.write("friction coefficient; \t \t \t " + str(self.a) +";	\t \t \t alpha; \n")
        f.write("simulation end time; \t \t \t " + str(90) +";	\t \t \t endTime; \n")
        f.close()       

    def Run(self):
        self.CreateDir()
        self.pro = sub.Popen([self.directory + "\\ODIS.exe"], shell=True)
        self.start = True

    def Wait(self):
        exit_code = self.pro.wait()
        self.complete = True
        return exit_code
        

    def IsRun(self):
        if self.start:
            if self.pro.poll() == None:
                return False
            else:
                self.complete = True
                return True
        else:
            return False

    def IsChild(self, proc, n):
         for pos in self.childrenPos:
            if pos == proc.node:
                return True
        
         return False

class Grid:
    def __init__(self, dim):
        self.dim = dim
        self.grid = [[0 for j in range(dim[0])] for i in range(dim[1])]
        self.max_proc = dim[0]*dim[1]

        self.running = []
        self.complete = []
        self.queue = []

    def PopulateGrid(self, hrange, arange):
        height = np.linspace(hrange[0],hrange[1],self.dim[1])
        alpha = np.linspace(arange[0],arange[1],self.dim[0])
        count = hcount = acount = 0
        for h in height:
            for a in alpha:
                self.grid[hcount][acount] = Process(count+1,h,a,(acount,hcount),self.dim)
                acount += 1
                count += 1

            acount = 0
            hcount += 1

    def RunProcess(self, p):
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
                                        self.RunProcess(newp)
          
                    self.CheckRunning()
                        
                else:
                    self.CheckRunning()

                looped+=1
                time.sleep(15)



            print("Queued processes: ", self.queue)
            print("Running processes: ", self.running)
            print("Completed processes: ", self.complete, "\n")

            
g = Grid((5,5))
pprint.pprint(g.grid)
g.PopulateGrid([100, 10000], [2.28e-7, 2.28e-8])
pprint.pprint(g.grid)
g.SolveGrid(3)


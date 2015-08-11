#! python

import subprocess as sub
import os
import shutil
import sys

class Process:
    def __init__(self, ID, height, alpha, Node, dim):
        self.id = ID
        self.h = height
        self.a = alpha
        self.start = False
        self.complete = False
        self.node = Node
        self.childrenPos = []
        self.directory = os.getcwd() + "\\h" + str(self.h) + "_alpha" + str(self.a)
        self.results_dir = self.directory #os.getcwd() + "\\Results\\h" + str(self.h) + "_alpha" + str(self.a)
        self.init = "true"

        if self.h <= 35:
            if self.node[0] == 0 and self.node[1] == dim[1] - 1: #top and right
                self.childrenPos.append((self.node[0]+1,self.node[1]))
                self.max_children = 1
            elif self.node[0] == 0: #right
                self.childrenPos.append((self.node[0]+1,self.node[1]))
                self.childrenPos.append((self.node[0],self.node[1]+1))
                self.max_children = 2
            elif self.node[0] == dim[0] - 1: #left
                self.max_children = 0
            else: #interior
                self.childrenPos.append((self.node[0]+1,self.node[1]))
                self.max_children = 1
        else:
            if self.node[0] == dim[0] - 1 and self.node[1] == 0: #bottom and left
                self.childrenPos.append((self.node[0],self.node[1]+1))
                self.max_children = 1
            elif self.node[1] == 0 and self.node[0] != dim[0] - 1: #bottom
                self.childrenPos.append((self.node[0]+1,self.node[1]))
                self.childrenPos.append((self.node[0],self.node[1]+1))
                self.max_children = 2
            elif self.node[1] == dim[1] - 1: #top
                self.max_children = 0
            else:
                self.childrenPos.append((self.node[0],self.node[1]+1))
                self.max_children = 1

        self.latnum = 90
        self.lonnum = 45
        self.dt = 100

        if (self.h < 35):
            self.latnum = 90 #2
        elif (self.h < 100):
            self.latnum = 60 # 3
        elif (self.h < 1000):
            self.latnum = 60 #3
        elif (self.h < 10000):
            self.latnum = 90 # 2
        elif (self.h <= 50000):
            self.latnum = 90 # 2
        else:
            self.latnum = 90 # 2
            self.dt = 25

        if (self.a >= 6e-8):
            self.dt = 105
        #elif (self.a >= 1e-8):
        #    self.dt = 64
        else:
            self.dt = 40
        #if (self.a < 1e-7):
        #    self.dt = 150
        #if (self.a < 5e-8):
        #    self.dt = 120
        #if (self.a < 3e-8):
        #    self.dt = 80
        #if (self.a <= 1e-8):
        #    self.dt = 50

        #if (self.h > 2800):
        #    self.dt = 120
        #if (self.h > 5000):
        #    self.dt = 100
        #if (self.h > 6500):
        #    self.dt = 80
        #''if (self.h > 8000):
         #   self.dt = 60
        #elif (self.h > 2800 and self.a < 5e-8):
        #    self.dt = 100
        #elif (self.h < 2800 and self.a < 3e-8):
        #    self.dt = 72

    def __repr__(self):
        return str(self.id)

    def __str__(self):
        return str(self.id)

    def CreateDir(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        # Always copy latest version of ODIS
        shutil.copy("ODIS.exe", self.directory)
        self.CreateInputFile()

    def CreateInputFile(self):
        f = open(self.directory+"\\input.in",'w')
        f.write("ocean thickness; \t \t \t " + str(self.h) +"; \t \t \t h; \n")
        f.write("friction coefficient; \t \t \t " + str(self.a) +"; \t \t \t alpha; \n")
        f.write("potential; \t \t \t OBLIQ; \t \t \t potential; \n")
        f.write("simulation end time; \t \t \t " + str(100) +"; \t \t \t endTime; \n")


        f.write("time step; \t \t \t " + str(self.dt) +"; \t \t \t timeStep; \n")
        f.write("latitude spacing; \t \t \t " + str(self.latnum) +"; \t \t \t lat; \n")
        f.write("longitude spacing; \t \t \t " + str(self.lonnum) +"; \t \t \t lon; \n")

        if self.node == (0,0):
            self.init = "false"
        f.write("initial conditions; \t \t \t " + self.init +"; \t \t \t init; \n")

        # elif (self.a > 1e-7):
        #     conv = 8e-7
        # elif (self.a < 1e-7 and self.h < 50):
        #     conv = 1e-6
        # elif (self.a > 1e-8 and self.h > 1000):
        #     conv = 5e-7
        # elif (self.a > 1e-8 and self.h < 1000):
        #     conv = 2e-7
        if (self.a <1e-7):
            conv = 5e-8
        elif (self.h < 100):
            conv = 5e-7
        else:
            conv = 1e-7
        #conv = 1e-7
        f.write("converge; \t \t \t " + str(conv) +"; \t \t \t convergence criteria; \n")


        f.close()

    def Run(self):
        self.CreateDir()
        if self.node == (0,0):
            input("Press any key after adding initial conditions to process 1...")
        print("attempting to run " + self.directory + "\\ODIS.exe")
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

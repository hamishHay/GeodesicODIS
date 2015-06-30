#! python

import subprocess as sub
import os
import shutil

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
        self.results_dir = os.getcwd() + "\\Results\\h" + str(self.h) + "_alpha" + str(self.a)

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
        return str(self.id)

    def __str__(self):
        return str(self.id)

    def CreateDir(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        shutil.copy("ODIS.exe", self.directory)
        self.CreateInputFile()

    def CreateInputFile(self):
        self.init = "true"
        f = open(self.directory+"\\input.in",'w')
        f.write("ocean thickness; \t \t \t " + str(self.h) +"; \t \t \t h; \n")
        f.write("friction coefficient; \t \t \t " + str(self.a) +"; \t \t \t alpha; \n")
        f.write("simulation end time; \t \t \t " + str(150) +"; \t \t \t endTime; \n")
        if self.node == (0,0):
            self.init = "false"
        f.write("initial conditions; \t \t \t " + self.init +"; \t \t \t init; \n")

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

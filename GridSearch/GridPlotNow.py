#! python
import Grid
import numpy as np

#a,h
g = Grid.Grid((10,80)) #10, 40
#h,a np.log10(50)
g.PopulateGrid([0,  4],[-5,-9]) # 2.69897
#g.SolveGrid(7)
g.CollectResults()

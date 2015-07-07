#! python
import Grid
import numpy as np

#a,h
g = Grid.Grid((8,50)) #10, 40
#h,a
g.PopulateGrid([np.log10(2), 3],[np.log10(5e-7),-9]) # 2.69897
g.SolveGrid(7)
g.CollectResults()

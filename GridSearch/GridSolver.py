#! python
import Grid
import numpy as np

#a,h
g = Grid.Grid((13,101)) #10, 40
#h,a np.log10(50)
g.PopulateGrid([0,  5],[-5,-9]) # 2.69897
#g.ResetRange(1e-10,1e-4,100,10001)
#g.MakeComplete(2e-9,8e-9,4000,10000)
g.SolveGrid(7,refine=False)
#g.CollectResults()
#g.PlotResults()

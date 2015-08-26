#! python
import Grid
import numpy as np

#a,h
g = Grid.Grid((25,71)) #10, 40
#h,a np.log10(50)
g.PopulateGrid([0,  4],[-5,-9]) # 2.69897
#####g.ResetRange(1e-10,1.1e-6,100,100000)
#g.MakeComplete(1.1e-7,1e-4,100,100000)
#g.SolveGrid(7,refine=False)
#g.DeleteRange(1e-10,2e-5,34.91,10001)
g.CollectResults()
g.PlotResults()

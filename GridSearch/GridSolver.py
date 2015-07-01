#! python
import Grid

g = Grid.Grid((10,40)) #10, 40
g.PopulateGrid([1.69897, 3], [1, 2.69897]) # 2.69897
g.SolveGrid(6)
g.CollectResults()

#! python
import Grid

g = Grid.Grid((8,40))
g.PopulateGrid([2, 3], [1, 2.69897]) # 2.69897
g.SolveGrid(6)
g.CollectResults()

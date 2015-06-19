#! python

#This is the master run file for exploration of h/alpha parameter space
#in an ODIS grid.



# 1. Receive input grid space to search over
# 2. Check current directory for nodes from previous runs that do not
#    need to be rerun.
# 3. Generate input files and new directories for starting nodes.
# 4. Execute user specified maximum number of nodes.
# 5. Routinely check each node for completion.
# 6. Copy output from completed node into new directory with new input files.
# 7. Begin 4 again until all nodes are complete.

import node
import inputs

init_data = inputs.Input()
init_data.ReadInputs("inputs.in")

nodes = []
nodes.append(node.Node(400,1e-7,init_data,False)




from ete3 import Tree
import sys

working_dir="/home/uvi/be/avs/store/dNdS_clones/"

program_name = sys.argv[0]
tree_name = sys.argv[1]
tree_dir = sys.argv[2]

if len(sys.argv) != 3:
    print ("Usage: prune_trees.py <input_tree> <tree_dir>")
    #sys.exit(1)

rt = Tree(tree_dir + "/all_clones/" + tree_name + ".raxml.bestTree")

#Prune tree maintaining clone leaves
#Example: clone_leaves = ["clone1", "clone2", "clone3"]

rt.prune(rt.get_leaves()[:-1])

rt.write(format=0, outfile =tree_dir + "/pruned/" + tree_name + "_pruned_rooted.nwk")

#Unroot the pruned tree
rt.unroot()
rt.write(format=0, outfile =tree_dir + "/pruned/" + tree_name + "_pruned_unrooted.nwk")
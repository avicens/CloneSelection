from ete3 import Tree
import sys

working_dir="/home/uvi/be/avs/store/dNdS_clones/prueba_COAD/"
tree_dir = working_dir + "trace/raxml_trees/"

program_name = sys.argv[0]
sample_name = sys.argv[1]
#Example: sample_name="TCGA-4N-A93T"

if len(sys.argv) != 2:
    print ("Usage: prune_trees.py <sample_name>")
    sys.exit(1)

#rt = Tree("/Users/avicens/Dropbox/proyectos/dNdS_clones/prueba_COAD/trace/raxml_trees/rooted/TCGA-4N-A93T_rooted.raxml.bestTree", format=0)
rt = Tree(tree_dir + "rooted/" + sample_name + "_rooted.raxml.bestTree")

#Prune tree maintaining clone leaves
#Example: clone_leaves = ["clone1", "clone2", "clone3"]

rt.prune(rt.get_leaves()[:-1])

#rt.write(format=0, outfile ="/Users/avicens/Dropbox/proyectos/dNdS_clones/prueba_COAD/trace/raxml_trees/pruned/TCGA-4N-93T_prunned.nwk")
rt.write(format=0, outfile =tree_dir + "pruned/" + sample_name + "_pruned_rooted.nwk")


#Unroot the pruned tree
rt.unroot()
rt.write(format=0, outfile =tree_dir + "pruned/" + sample_name + "_pruned_unrooted.nwk")
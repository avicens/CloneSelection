from ete3 import Tree

rt = Tree("/Users/avicens/Dropbox/proyectos/dNdS_clones/prueba_COAD/trace/raxml_trees/rooted/TCGA-4N-A93T_rooted.raxml.bestTree", format=0)

#Prune tree maintaining clone leaves

#Ej: clone_leaves = ["clone1", "clone2", "clone3"]

rt.prune(rt.get_leaves()[:-1])

rt.write(format=0, outfile ="/Users/avicens/Dropbox/proyectos/dNdS_clones/prueba_COAD/trace/raxml_trees/pruned/TCGA-4N-93T_prunned.nwk")

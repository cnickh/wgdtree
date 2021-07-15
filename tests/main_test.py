from wgdtree import sim, root
from ete3 import Tree, PhyloTree

NUM_TREES = 100
trees = ["balanced.newick","caterpillar.newick","small.newick"]

for tree in trees:
    s_tree = PhyloTree("../example_trees/"+tree,format=1)

    sim(s_tree,NUM_TREES,tree+"_trees")


    for num in range(NUM_TREES):

        g_tree = PhyloTree(tree+"_trees/tree"+ str(num) +".newick",format=1)

        g_tree = root(g_tree,s_tree)



        filename = "finished_"+tree+"/tree" + str(num) + ".newick.labeled"
        g_tree.write(format=0,outfile=filename,features=["evoltype"])


        #g_tree = wgdtree.rrate.find_WGD(s_tree,g_tree)

        #filename = "finished_"+tree+"/tree" + str(num) + ".newick.labeled"
        #g_tree.write(format=0,outfile=filename,features=["event0","event1"])

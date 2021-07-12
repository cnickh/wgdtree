import sys
from ete3 import PhyloTree
#from wgdtree.split_it import split_it, break_up
from wgdtree.root import reconcile

def sortFunc(e):
    return e['len']

#reduce tree size for processing 
def break_up(tree,size=100):
    
    nw_trees =[]
    poss_trees = []
    
    while(len(tree.get_leaves()) > size):  
    
        poss_trees.clear()
        for n in tree.traverse():
            if len(n.get_leaves()) < size:
                poss_trees.append({'node': n, 'len': len(n.get_leaves())})
                
        poss_trees.sort(key=sortFunc)
        _n = poss_trees.pop()['node']
        print("Removing: " + str(len(_n.get_leaves())) + " | Tree_total: " + str(len(tree.get_leaves())))
        nw_trees.append(_n.detach())
        
        
        
    nw_trees.append(tree)
    return nw_trees
    


    

gene_trees = []
gene_trees.append(PhyloTree("labelednalignment14.nwk",format=1))
gene_trees.append(PhyloTree("labelednalignment20.nwk",format=1))
s_tree = PhyloTree("moderntry.nwk",format=1)

for t in gene_trees:
    t.set_outgroup(t.children[0])
    reconcile(t,s_tree)
    print(len(break_up(t)))


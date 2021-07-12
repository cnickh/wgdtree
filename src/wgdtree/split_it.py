from ete3 import Tree


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
    

    
    

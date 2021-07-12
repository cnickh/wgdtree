from ete3 import Tree
from ete3 import PhyloTree
import numpy as np



def GET_SPECIES(n):
    return n.name.split('_')[0]



def find_dups(gene_tree):

    children = []
    dups = 0
       
    for node in gene_tree.traverse():
           
        children = node.get_children()
        
        if children:
            for n in children[0].get_leaves():
                species1 = GET_SPECIES(n)
                
                for q in children[1].get_leaves():
                    species2 = GET_SPECIES(q)
                        
                    if(species1 == species2):
                        node.add_feature('evoltype','D')
                        dups += 1
                        break
            
                if 'evoltype' in node.features: 
                    break
            
            if('evoltype' not in node.features):
                node.add_feature('evoltype','S')
            
            
            
    return dups



def lca(G,S): # map gene tree node to species tree
    
    species = []
    s_match = []

    #get species present in gene tree 'G'
    for n in G.get_leaves():
        s = GET_SPECIES(n)
                    
        if s not in species:
            species.append(s)
            
            #match species on 'G' to a node on  'S'
            for l in S.get_leaves():
                if s in l.name:
                    s_match.append(l)

    if len(s_match) > 1:

        try:
            return S.get_common_ancestor(s_match)
        except:
            #print(s_match)
            #print(S)
            #exit()
            return S

    else:
            
        return s_match[0]




def reconcile(G,_S):
    
    S = _S.copy()
    s_children = S.get_children()
    g_children = G.get_children()
    
    if not g_children and not s_children:
        return G
    
    #1
    if s_children:
        if(lca(G,S) in s_children[0].traverse()):
            s_children[0].up.add_child(reconcile(G,s_children[0]))
            s_children[0].detach() #replace
            return S
            
        elif(lca(G,S) in s_children[1].traverse()):
            s_children[1].up.add_child(reconcile(G,s_children[1]))
            s_children[1].detach() #replace
            return S
    
    #2
    if g_children and s_children:
        if(lca(G,S) == S 
           and (lca(g_children[0],S) in s_children[0].traverse()) 
           and (lca(g_children[1],S) in s_children[1].traverse())
           ):
            t = PhyloTree() #join
            t.add_child(reconcile(g_children[0],s_children[0]))
            t.add_child(reconcile(g_children[1],s_children[1]))
            return t
            
        elif(lca(G,S) == S 
             and lca(g_children[0],S) in s_children[1].traverse() 
             and lca(g_children[1],S) in s_children[0].traverse()
             ):
            t = PhyloTree() #join
            t.add_child(reconcile(g_children[0],s_children[1]))
            t.add_child(reconcile(g_children[1],s_children[0]))
            return t

    #3
    if g_children:
        if(lca(G,S)==S 
           and (lca(g_children[0],S) == S  or lca(g_children[1],S) == S)
           ):
            t = PhyloTree() #join
            t.add_child(reconcile(g_children[0],S))
            t.add_child(reconcile(g_children[1],S))
            return t
    
    #print("Houston we have a problem")
    #print("G:")
    #print(G)
    #print(lca(G,S))
    #if g_children:
        #print(lca(g_children[0],S))
        #print(lca(g_children[1],S))
    #print("S:")
    #print(S)
    print("killed on err")
    exit()
    
    
def root(gene_tree,species_tree):
        
    minimun = 999
        
    gene_tree.unroot()
            
    for node in gene_tree.children():
        
        if node != gene_tree:
            
            tree_poss = gene_tree.copy()
            
            try:
                tree_poss.set_outgroup(tree_poss&node.name)
            except:
                print("Assert ERR: outgroup not set " + node.name)
                
            if(tree_poss): 
                
                dups = find_dups(reconcile(tree_poss,species_tree))
                
                if(dups < minimun):
                    best_root = node
                    minimun = dups

    #try:
    gene_tree.set_outgroup(best_root)
    tree = reconcile(gene_tree,species_tree)
    find_dups(tree)
    return tree
    #except:
        #print("error!!")   
        #return


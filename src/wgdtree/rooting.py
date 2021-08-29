from ete3 import Tree
from ete3 import PhyloTree
import numpy as np



def GET_SPECIES(n):
    return n.name.split('_')[0]

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
            return S

    else:
            
        return s_match[0]


def find_loss(gene_tree):
    losses = 0
    
    for l in gene_tree.get_leaves():
        if "_" not in l.name:
            l.name = l.name + "_LOST"
            losses+=1
    
    return losses


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

def join(n1,n2):
    
    n0 = PhyloTree()
    n0.name = "in_node"
    
    if n1.dist == n2.dist:
        n0.dist = n1.dist/2
        n1.dist = n1.dist/2
        n2.dist = n2.dist/2
    

    else:
        n0.dist = 25
    
    n0.add_child(n1)
    n0.add_child(n2)
            
    return n0


def reconcile(G,_S):
    
    S = _S.copy()
    s_children = S.get_children()
    g_children = G.get_children()
    
    if not g_children and not s_children:
        return G
    
    #1
    if s_children:
        if(lca(G,S) in s_children[0].traverse()):
            parent = s_children[0].up
            s_children[0].detach() #replace
            parent.add_child(reconcile(G,s_children[0]))
            return S
            
        elif(lca(G,S) in s_children[1].traverse()):
            parent = s_children[1].up
            s_children[1].detach() #replace
            parent.add_child(reconcile(G,s_children[1]))
            return S
    
    #2
    if g_children and s_children:
        if(lca(G,S) == S 
           and (lca(g_children[0],S) in s_children[0].traverse()) 
           and (lca(g_children[1],S) in s_children[1].traverse())
           ):

            n1 = reconcile(g_children[0],s_children[0])
            n2 = reconcile(g_children[1],s_children[1])
            return join(n1,n2)
            
        elif(lca(G,S) == S 
             and lca(g_children[0],S) in s_children[1].traverse() 
             and lca(g_children[1],S) in s_children[0].traverse()
             ):

            n1 = reconcile(g_children[0],s_children[1])
            n2 = reconcile(g_children[1],s_children[0])
            return join(n1,n2)

    #3
    if g_children:
        if(lca(G,S)==S 
           and (lca(g_children[0],S) == S  or lca(g_children[1],S) == S)
           ):
            
            n1 = reconcile(g_children[0],S)
            n2 = reconcile(g_children[1],S)
            return join(n1,n2)
    
    print("Reconciliation failed! Maybe your tree isn't binary")
    exit()

    
def root(gene_tree,species_tree):
        
    dups = 999999

    for root in gene_tree.get_children():
        tree = gene_tree.copy()
        tree.set_outgroup(tree&root.name)
        r_tree = reconcile(tree,species_tree)
        temp = find_dups(r_tree) + find_loss(r_tree)
        if(temp < dups):
            best_root = root
            dups = temp

    gene_tree.set_outgroup(best_root)
    tree = reconcile(gene_tree,species_tree)
    find_dups(tree)
    find_loss(tree)
    return tree


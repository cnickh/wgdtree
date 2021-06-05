
from ete3 import Tree, TreeStyle

#places wgd events on the gene tree
def find_WGD(species_root,gene_tree,gene_wgd,wgd_num):

    #get children
    children = species_root.get_children()
                        
    #check if were at a leaf
    if children:
        nodeR = children[0]
        nodeL = children[1]
        
        #check if node is a duplication
        if('WGD' in species_root.features and species_root.WGD == 'Y'): 
            pSpecies = []
                
            #get all species names present for wgd
            for p in species_root.get_leaves():
                pSpecies.append(p.name)
                    
            if (pSpecies):    
                cmnleaves = [] #create array to match trees
                
                #search gene_tree for all leaves that corespond to affected species
                for g in gene_wgd.get_leaves():
                    for s in pSpecies:
                        if s in g.name:
                            cmnleaves.append(g.name)
                                    
                if(cmnleaves):
                             
                    #create event name         
                    eventname = "event" + str(wgd_num)  
                             
                    #find matching node in gene tree, .get_mrca_idx_from_tip_labels() 
                    gene_wgd = gene_tree.idx_dict[gene_tree.get_mrca_idx_from_tip_labels(cmnleaves)]

                    #place/label duplication
                    if('D' in gene_wgd.features and gene_wgd.D == 'Y'):
                        gene_wgd.add_feature(eventname,"P")
                        
                        children = gene_wgd.get_children()
                        
                        if(children):   
                        
                            find_WGD(nodeR,gene_tree,children[0],wgd_num+1)
                            find_WGD(nodeL,gene_tree,children[1],wgd_num+1)
                        
                            find_WGD(nodeL,gene_tree,children[0],wgd_num+1)
                            find_WGD(nodeR,gene_tree,children[1],wgd_num+1)
            
                        return
                        
                    else:
                        gene_wgd.add_feature(eventname,"M")
                        wgd_num+=1
                                             
        find_WGD(nodeR,gene_tree,gene_wgd,wgd_num)
        find_WGD(nodeL,gene_tree,gene_wgd,wgd_num)
    
    return




#count losses 
def rrates(gene_tree,pSpecies):
    
    losses = 0
    #find most recent wgd and get all present sequences
    for g in gene_tree.treenode.traverse():
        if("event0" in g.features):   
            for n in g.traverse():
                if ("event1" in n.features):
                    copies = 0
                    for l in n.get_leaves():
                        for s in pSpecies:
                            if s in l.name and "*LOST" not in l.name:
                                copies+=1
                    if copies > 2*len(pSpecies):
                        losses += 0
                    else:
                        losses += 2*len(pSpecies)-copies
                                
            return (losses,g.event0)
                                                                        
    return 0

#Actually doing stuff starts here#

num=0
lpossible = 0
rpossible = 0
total_llosses = 0
total_rlosses = 0

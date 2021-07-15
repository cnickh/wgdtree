from ete3 import Tree


#@desc this function will make an attempt to attribute duplications labelled on a gene tree to
# wgd events on a species tree.
def place_wgd(species_root,gene_wgd,wgd_num=0):

    #get children        
        
    #check if node is a duplication
    if('WGD' in species_root.features and species_root.WGD == 'Y'): 
        pSpecies = []
            
        #get all species names present for wgd
        if len(species_root.get_children()) < 2:
            pSpecies.append(species_root.name)
        else:
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
                             
                #find matching node in gene tree
                try:
                    gene_wgd = gene_wgd.get_common_ancestor(cmnleaves)
                except:
                    gene_wgd = gene_wgd

                    
                #check if were at the correct node make recursive call if were at the wrong node
                for node in gene_wgd.get_leaves():
                    if node.name.split('_')[0] not in pSpecies:
                        children = gene_wgd.get_children()
                        if(children):
                                
                            place_wgd(species_root,children[0],wgd_num)
                            place_wgd(species_root,children[1],wgd_num)
                                
                        return
                
                children = gene_wgd.get_children()
                if children:
                    num = 0
                    for c in children:
                        if 'evoltype' in c.features and c.evoltype == 'D':
                            num +=1
                            
                    if num == 2:
                        place_wgd(species_root,children[0],wgd_num)
                        place_wgd(species_root,children[1],wgd_num)
                        return
                    

                #place/label duplication
                if('evoltype' in gene_wgd.features and gene_wgd.evoltype == 'D'):
                    gene_wgd.add_feature(eventname,"P")
                        
                    children = species_root.get_children()
                    if children:
                        nodeR = children[0]
                        nodeL = children[1]
                        
                        children = gene_wgd.get_children()
                            
                        if(children):   
                            
                            place_wgd(nodeR,children[0],wgd_num+1)
                            place_wgd(nodeL,children[1],wgd_num+1)
                            
                            place_wgd(nodeL,children[0],wgd_num+1)
                            place_wgd(nodeR,children[1],wgd_num+1)
            
                    return
                        
                else:
                    gene_wgd.add_feature(eventname,"M")
                    wgd_num+=1
                    
    children = species_root.get_children()
    if children:
        nodeR = children[0]
        nodeL = children[1]
        place_wgd(nodeR,gene_wgd,wgd_num)
        place_wgd(nodeL,gene_wgd,wgd_num)                                        

    return 0




#count copies 
def rrates(gene_tree,pSpecies):
    
    p=0
    for l in gene_tree.get_leaves():
        if pSpecies in l.name and "*LOST" not in l.name:
            p+=1
    if p==0:
        return 0
    
    copies = 0
    #find most recent wgd and get all present sequences
    for g in gene_tree.traverse():
        if("event0" in g.features):   
            for n in g.traverse():
                if ("event1" in n.features):
                    for l in n.get_leaves():
                        if pSpecies in l.name and "*LOST" not in l.name:
                            copies+=1
                                
            return (copies,g.event0)
                                                                        
    return -1

#Actually doing stuff starts here#
def rrates_batch(gene_trees,species_tree,pSpecies):
    num=0
    lpossible = 0
    rpossible = 0
    total_llosses = 0
    total_rlosses = 0
    for tree in gene_trees:
        
        place_wgd(species_tree,tree)
        
        filename = "finished_trees/tree" + str(num) + ".newick.labeled"
        tree.write(filename,format=0,features=["event0","event1"])
        
        results = rrates(tree,pSpecies)

        if results[1] == "M":
            lpossible += len(pSpecies)*2
            total_llosses += results[0]
            
        elif results[1] == "P":
            rpossible += len(pSpecies)*4
            total_rlosses += results[0]
        num+=1
    if(lpossible):
        print("Percent lossed after loss: " + str(total_llosses)+"/"+str(lpossible))
    print("Percent lossed after a retained event: " + str(total_rlosses)+"/"+str(rpossible))
    
    return

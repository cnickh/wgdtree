from ete3 import Tree

#@desc this function will make an attempt to attribute duplications labelled on a gene tree to
# wgd events on a species tree.
def place_wgd(species_root,gene_wgd,wgd_num=0):      
    
    children = species_root.get_children()
    
    if('WGD' in species_root.features and species_root.WGD == 'Y'):
                
        #match species tree node to gene tree
        species = []
        cmn_leaves = []
        
        for s in species_root.get_leaves():
            species.append(s.name)
            
            for g in gene_wgd.get_leaves():
                if s.name in g.name:
                    cmn_leaves.append(g)
        
        #if we were able to match leaves to the gene tree
        if cmn_leaves:
            
            nxt_wgd = gene_wgd.get_common_ancestor(cmn_leaves)
            
            leaves = []
            
            #get all the species under the chosen node on the gene tree
            for l in nxt_wgd.get_leaves():
                if l.name.split('_')[0] not in leaves:
                    leaves.append(l.name.split('_')[0])
            
            #make sure were in the right place before placing

            if sorted(species) == sorted(leaves): 
                
                c = nxt_wgd.get_children()
                
                #place event
                eventname = "event" + str(wgd_num)
                
                if 'evoltype' in nxt_wgd.features and nxt_wgd.evoltype == 'D':
                       
                    #make sure we are placing wgd events to minimize ssd events
                    if c and 'evoltype' in c[0].features and 'evoltype' in c[1].features and c[0].evoltype == 'D' and c[1].evoltype == 'D':
                        
                        place_wgd(species_root,c[0],wgd_num)
                        place_wgd(species_root,c[1],wgd_num)
                
                    else:
                        
                        nxt_wgd.add_feature(eventname,"P")
                            
                        #make recursive calls
                        if children and c:
                                                            
                            place_wgd(children[0],c[0],wgd_num+1)
                            place_wgd(children[1],c[1],wgd_num+1)
                                                            
                            place_wgd(children[1],c[0],wgd_num+1)
                            place_wgd(children[0],c[1],wgd_num+1)

                else:
                    for l in nxt_wgd.get_leaves():
                        if "_LOST" not in l.name:
                            nxt_wgd.add_feature(eventname,"M")
                            break
                    
                    if children:
                        place_wgd(children[0],gene_wgd,wgd_num+1)
                        place_wgd(children[1],gene_wgd,wgd_num+1)
                
            elif len(leaves) > len(species):
                place_wgd(species_root,gene_wgd.children[0],wgd_num)
                place_wgd(species_root,gene_wgd.children[1],wgd_num)
                
    else:
        if children:
            place_wgd(children[0],gene_wgd,wgd_num)
            place_wgd(children[1],gene_wgd,wgd_num)
    
    return 0

#count copies 
def rrates(gene_tree,pSpecies):
    
    #initialize variables
    rcopies = 0
    lcopies = 0
    l_poss = 0
    r_poss = 0
    p = 0
    
    for l in gene_tree.get_leaves():
        if pSpecies in l.name and "_LOST" not in l.name:
            p+=1

    if p!=0:
         
        for g in gene_tree.iter_descendants():
            if("event0" in g.features):
                copies = 0
                c_poss = 0
                
                for n in g.iter_descendants():
                    if "event1" in n.features:
                        c_poss += 2
                        temp = 0
                        for l in n.get_leaves():
                            if pSpecies in l.name and "_LOST" not in l.name:
                                temp+=1
                                
                        if temp > 2:
                            copies += 2
                        else:
                            copies += temp
                
                if g.event0 == 'M':
                    lcopies += copies
                    l_poss += c_poss
                
                elif g.event0 == 'P':
                    rcopies += copies
                    r_poss += c_poss

    if r_poss == 0 and l_poss == 0:
        ret = -1
    elif rcopies == 0 and lcopies == 0 and p != 0:
        ret = -1
    else:
        ret = ((rcopies,r_poss),(lcopies,l_poss),p)
    
    return ret

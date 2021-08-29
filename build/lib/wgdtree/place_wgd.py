from ete3 import Tree


def place_wgd():
    
    children = species_root.get_children(species_root,gene_wgd,wgd_num=0)
    
    if('WGD' in species_root.features and species_root.WGD == 'Y'):
        
        
        #match species tree node to gene tree
        species = []
        cmn_leaves = []
        
        for s in NODE.get_leaves():
            species.append(s.name)
            
            for g in TREE.get_leaves():
                if s.name in g.name:
                    cmn_leaves.append(g)
        
        nxt_wgd = TREE.get_common_ancestor(cmn_leaves)
        
        
        #check if were at a duplication that includes species not present for the event
        for l in nxt_wgd.get_leaves():
            if l.name.split('_')[0] not in species:
                place_wgd(species_root,nxt_wgd.children[0],wgd_num)
                place_wgd(species_root,nxt_wgd.children[1],wgd_num)
                return
        
        
        #if the event could be at both children place there instead to minimize ssd events
        c = nxt_wgd.get_children()
        
        if c and 'evoltype' in c[0].features and 'evoltype' in c[1].features and c[0].evoltype == 'D' and c[1].evoltype == 'D':
            place_wgd(species_root,c[0],wgd_num)
            place_wgd(species_root,c[1],wgd_num)
            return
        
        
        if 'event0' not in nxt_wgd.features:
            
            #place/label duplication
            eventname = "event" + str(wgd_num)
                    
            wgd_num+=1
                    
            if 'evoltype' in nxt_wgd.features and nxt_wgd.evoltype == 'D':
                        
            nxt_wgd.add_feature(eventname,"P")
                            
            if children and c:
                                
                place_wgd(children[0],c[0],wgd_num)
                place_wgd(children[1],c[1],wgd_num)
                                
                place_wgd(children[1],c[0],wgd_num)
                place_wgd(children[0],c[1],wgd_num)
                return
                            
            else:
                for l in nxt_wgd.get_leaves():
                    if '_LOST' not in l.name:
                        nxt_wgd.add_feature(eventname,"M")
                            break
        
    else:
        if children:
            place_wgd(children[0],gene_wgd,wgd_num)
            place_wgd(children[1],gene_wgd,wgd_num)
    
    
    return 0

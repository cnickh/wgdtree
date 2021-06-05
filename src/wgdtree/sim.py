import sys
from ete3 import Tree, TreeStyle
from random import randint
from random import uniform
from scipy.stats import expon


#=================================Helper=Functions==============================
#function for creating duplication in a tree
def dup_root(g_root,t):
    
        
    cpy = g_root.copy("newick-extended")
    cpy0 = g_root.copy("newick-extended")
    
    g_root.name = g_root.name + "_"    
        
    #fix branch dist issue
    (cpy).dist = g_root.dist-t
    (cpy0).dist = g_root.dist-t
    g_root.dist = t
            
    #ensure only one wgd feature 
    if("WGD" in cpy.features):
        (cpy).del_feature("WGD")
        (cpy0).del_feature("WGD")

    
    children = g_root.get_children()
    
    if children:
        #remove previous children
        children[0].detach()
        children[1].detach()
            
    #add duplicates
    g_root.add_child(cpy)
    g_root.add_child(cpy0)
        
    return

#function for creating losses in a tree
def lose_root(g_root):
        
    
    while len(g_root.up.get_children()) == 1:
        if(g_root.up):
            g_root = g_root.up
                
    g_root.detach()
            
    return

#this allows trees to be read by notung
def remove_single_links(root):
    
    children = root.get_children()
    
    if len(children) == 1:
        children[0].dist += root.dist
        parent = root.up
        parent.add_child(children[0])
        remove_single_links(children[0])
        root.detach()
        
    elif len(children) == 2:
        remove_single_links(children[0])
        remove_single_links(children[1])
    
    return
#===============================================================================

#Simulation
def bd_model(root,g_tree,t,c):
        
    global loss_rate
    global ssd_rate
            
    wgd_t = 999

    if c > 1:
        rate_sum = loss_rate + ssd_rate
    else:
        rate_sum = ssd_rate

    #create specified wgd events
    if("WGD" in root.features and root.WGD == "Y"):
        wgd_t = int(root.T)

    #get time for first event
    nxt_event = expon.rvs()/(rate_sum*c)
    branch = root.dist - nxt_event
    
    #check for wgd
    if wgd_t < nxt_event:
        wgd_t = 999
        nxt_event -= wgd_t
        c=c*2
        for n in g_tree.search_nodes(name=root.name):
            dup_root(n,wgd_t)
            n.add_feature("WGD","Y")
    
    while(branch > 0):
    
        #set rate_sum
        if c > 1:
            rate_sum = loss_rate + ssd_rate
        else:
            rate_sum = ssd_rate
            
        #time for some loss_rate:ssd_rate sorcery
        draw0 = uniform(0,1)
        if(draw0 < ssd_rate/rate_sum):
            c+=1
            #pick random copy to duplicate
            g_poss = g_tree.search_nodes(name=root.name)
            g_root = g_poss[randint(0,len(g_poss)-1)] 
            #duplicate
            dup_root(g_root,nxt_event)
        else:
            c-=1
            #pick random copy to lose
            g_poss = g_tree.search_nodes(name=root.name)
            g_root = g_poss[randint(0,len(g_poss)-1)]            
            #lose
            lose_root(g_root)
            
        nxt_event = expon.rvs()/(rate_sum*c)
        branch -= nxt_event
        
        #do some magic to add wgd
        if wgd_t < root.dist - branch:
            wgd_t = 999
            time_of_event = (root.dist - branch) - wgd_t
            nxt_event -= time_of_event
            branch += time_of_event
            c=c*2
            for n in g_tree.search_nodes(name=root.name):
                dup_root(n,time_of_event)
                n.add_feature("WGD","Y")

    children = root.get_children()

    if children:
        bd_model(children[0],g_tree,t,c)
        bd_model(children[1],g_tree,t,c)

    return



def sim(species_tree,gene_tree,n):
    for int i in range(0,n):
        
    return







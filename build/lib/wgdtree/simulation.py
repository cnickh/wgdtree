from ete3 import Tree#, TreeStyle
from random import randint
from random import uniform
from scipy.stats import expon


global ssd_rate
global loss_rate

#=================================Helper=Functions==============================
#function for creating duplication in a tree
def dup_root(g_root,t):
    
        
    cpy = g_root.copy("newick-extended")
    cpy0 = g_root.copy("newick-extended")
    
    g_root.name = g_root.name + "_"    
        
    #fix branch dist issue
    cpy.dist = g_root.dist-t
    cpy0.dist = g_root.dist-t
    g_root.dist = t
            
    #ensure only one wgd feature 
    if("WGD" in cpy.features):
        cpy.del_feature("WGD")
        cpy0.del_feature("WGD")

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
def bd_model(root,g_tree,c):
        
    global loss_rate
    global ssd_rate
            
    wgd_t = 99999999

    if c > 1:
        rate_sum = loss_rate + ssd_rate
    else:
        rate_sum = ssd_rate

    #create specified wgd events
    if("WGD" in root.features and root.WGD == "Y"):
        wgd_t = int(root.T)
    
    branch = root.dist

    while(branch > 0):
    
        #get time of next event
        nxt_event = expon.rvs()/(rate_sum*c)

        #do some magic to add wgd
        if wgd_t < nxt_event + (root.dist-branch) and wgd_t >= (root.dist - branch):
            
            time_of_event = wgd_t - (root.dist - branch)
            branch -= time_of_event
            
            nodes_to_dup = g_tree.search_nodes(name=root.name)
            for n in nodes_to_dup:
                c+=1
                dup_root(n,time_of_event)
                n.add_feature("WGD","Y")
                
            wgd_t = 9999999


        else:
            if branch - nxt_event > 0:
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
                
            branch -= nxt_event

    children = root.get_children()

    if children:
        bd_model(children[0],g_tree,c)
        bd_model(children[1],g_tree,c)

    return

def sim(species_tree,n,path=".",ssd=.002,loss=.002):
    
    global ssd_rate
    global loss_rate
    
    
    ssd_rate = ssd
    
    loss_rate = loss
    
    
    trees = []
    for x in range(n):
        trees.append(species_tree.copy("newick-extended"))
        for q in trees[x].traverse():
            if("WGD" in  q.features):
                q.del_feature("WGD")

    num = 0

    for x in range(n):
                
        #build gene_trees and simulate losses 
        bd_model(species_tree,trees[x],1)
        remove_single_links(trees[x])

        #randomize names
        for n in trees[x].get_leaves():
            n.name = n.name + "_" + str(uniform(0,9999999))

        #write trees to files
        filename = path + "/tree" + str(num) + ".newick"
        trees[x].write(format=0,outfile=filename,features=["WGD"])
        num +=1
        
    return 0
    





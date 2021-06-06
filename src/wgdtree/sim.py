import sys
from ete3 import Tree#, TreeStyle
from random import randint
from random import uniform
from scipy.stats import expon

NUM_TREES = 100 

#readin labelled species tree
in_stree = sys.argv[1]
loss_rate = float(sys.argv[2])
ssd_rate = float(sys.argv[3])
sf = open(in_stree, "r")
species_tree = Tree(sf.read(),format=1)

#copy basic topology of species tree for gene trees
trees = []
for x in range(NUM_TREES+1):
    trees.append(species_tree.copy("newick-extended"))

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
            wgd_t = 9999999
            nxt_event -= time_of_event
            branch -= time_of_event
            for n in g_tree.search_nodes(name=root.name):
                c+=1
                dup_root(n,time_of_event)
                n.add_feature("WGD","Y")


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
                if not len(g_poss):
                    print("No g_poss")
                g_root = g_poss[randint(0,len(g_poss)-1)] 
                #duplicate
                dup_root(g_root,nxt_event)
            else:
                c-=1
                if c == 0:
                    print("nigga")
                #pick random copy to lose
                g_poss = g_tree.search_nodes(name=root.name)
                if not len(g_poss):
                    print("No g_poss")
                g_root = g_poss[randint(0,len(g_poss)-1)]            
                #lose
                lose_root(g_root)
            
        branch -= nxt_event

    children = root.get_children()

    if children:
        bd_model(children[0],g_tree,t,c)
        bd_model(children[1],g_tree,t,c)

    return


#call everything
num = 0

#print(randint(0,-1))

for x in range(NUM_TREES):
    
    #build gene_trees and simulate losses 
    bd_model(species_tree,trees[x],0,1)
    remove_single_links(trees[x])

    #randomize names
    for n in trees[x].get_leaves():
        n.name = n.name + "_" + str(randint(0,99999))

    #write trees to files
    filename = "raw_sim/tree" + str(num) + ".newick"
    trees[x].write(format=0,outfile=filename,features=["WGD"])
    num +=1


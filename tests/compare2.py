from ete3 import Tree
from os import listdir


NUM_TREES = 500

trees = ["balanced.newick","caterpillar.newick","small.newick"]

failed_trees = []

raw_trees = []
labeled_trees = []

for x in range(NUM_TREES):
    raw_trees.append(Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);"))
    labeled_trees.append(Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);"))




def compare(root0,root1,num):
    global total_events
    global correctly_placed
        
        
    for n in root0.traverse():
        
        if ('event0' in n.features and n.event0 == 'P') or ('event1' in n.features and n.event1 == 'P'):
            
            cmnleaves = []
            
            for g in n.get_leaves():
                cmnleaves.append(g.name)   
                
            n2 = root1.get_common_ancestor(cmnleaves)
            
            if 'WGD' in n2.features:
                correctly_placed += 1
                total_events += 1
            else:
                if num not in failed_trees:
                    failed_trees.append(num)
                total_events += 1
    return



for tree in trees:
    print(tree)
    for p in listdir(tree+"_trees/"):
            i = int(p.split('.')[0][4:])
            f = open(tree+"_trees/"+p, "r")
            raw_trees[i] = Tree(f.read(),format=0)
        
    for p in listdir("finished_"+tree+"/"):
        i = int(p.split('.')[0][4:])
        f = open("finished_"+tree+"/"+p, "r")
        labeled_trees[i] = Tree(f.read(),format=0)

    total_events = 0
    correctly_placed = 0
    failed_trees = []
    for i in range(NUM_TREES):
        compare(labeled_trees[i],raw_trees[i],i)
        
    print(failed_trees)
    print(str(correctly_placed) + " of " + str(total_events) + " placed correctly.")

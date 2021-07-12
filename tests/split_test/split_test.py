import sys
from ete3 import Tree
#from wgdtree.split_it import split_it, break_up

#define our recursive function
def split_it(root,tosplit=0):
        
    #get children
    children = root.get_children()
    
    nw_trees = []

    #check if were at a leaf
    if children:
        
        #check if node is a duplication
        if('evoltype' in root.features and root.evoltype == 'D'): 
                        
            for x in split_it(children[0],1):
                nw_trees.append(x)
            for x in split_it(children[1],1):
                nw_trees.append(x)            

            
        elif(tosplit):
            
            i=0
            for x in split_it(children[0],0):
                nw_trees.append(x)
                i+=1
        
            q=0  
            for x in split_it(children[1],0):
                nw_trees.append(x)
                q+=1
                
            if not i and q:
                children[0].detach()
                nw_trees.append(children[0])
                
            elif not q and i:
                children[1].detach()
                nw_trees.append(children[1])

            elif not nw_trees:
                root.detach()
                nw_trees.append(root)

    return nw_trees


#set filenames do set-up etc...
gene_batch = sys.argv[1]

#readin reconciled gene trees and species tree
gf = open(gene_batch, "r")

gene_trees = []
names =[]
line = gf.readline()
while(line): 
    f = open(line[:-1], "r")
    gene_trees.append(Tree(f.read(),format=1))
    line = gf.readline()
    names.append(line)


#call our function
i=1
for tree in gene_trees:
    print(split_it(tree))
    i+=1

#end

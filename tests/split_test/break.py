import sys
from ete3 import Tree


def sortFunc(e):
    return e['len']

#reduce tree size for processing 
def break_up(tree,size=100):
    
    nw_trees =[]
    poss_trees = []
    
    while(len(tree.get_leaves()) > size):  
    
        poss_trees.clear()
        for n in tree.traverse():
            if len(n.get_leaves()) < size:
                poss_trees.append({'node': n, 'len': len(n.get_leaves())})
                
        poss_trees.sort(key=sortFunc)
        _n = poss_trees.pop()['node']
        nw_trees.append(_n.detach())
        
        
        
    nw_trees.append(tree)
    return nw_trees
    


#set filenames do set-up etc...
gene_batch = sys.argv[1]

#readin trees
gf = open(gene_batch, "r")

gene_trees = []
names =[]
line = gf.readline()
while(line): 
    f = open(line[:-1], "r")
    gene_trees.append(Tree(f.read(),format=1))
    line = gf.readline()
    names.append(line[:-1])


#call our function
i = 0
for tree in gene_trees:
    
    sub_trees = break_up(tree)
    q=0
    for t in sub_trees:
        filename = "sub_trees/" + names[i]+".sub_tree" +str(q)
        t.write(format=0,outfile=filename,features=[])
        q+=1
        
    i +=1


#end

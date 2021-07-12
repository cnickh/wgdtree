# wgdtree Package

This package provides the following functions:


Function | Description
-------- | -----------
rrates(labeled_gene_tree) |  returns conditional probability of retention
find_wgd(species_tree, gene_tree) |  returns gene tree with labeled wgd events
sim(species_tree,n,path) | returns 0, writes 'n' simulated trees to 'path'
root(species_tree,gene_tree,mode) | roots gene_tree with species_tree, roots vary by mode        
        
*all trees must be in .newick format, and the ete3 package is required.

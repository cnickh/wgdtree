# wgdtree Package


install with ```$ pip install wgdtree ```

```python
import wgdtree
from ete3 import PhyloTree
````

This package provides the following functions:


**sim(species_tree,n,path,ssd_rate,loss_rate)** | This function accepts a species tree, an integer value 'n', a path, a loss rate floating point and a ssd rate also floating point. The function will simulate n statistically probable gene trees and then write them as newick files to the given path. This function allows for wgd duplications that are labelled on the species tree to be accounted for. Annotations are expected to be in NHX format.


```python
path="/home/user/simulated_trees_dir"

n = 10

ssd_rate = 0.002

loss_rate = 0.002

species_tree = PhyloTree("example_species_tree.newick")

wgdtree.sim(species_tree,n,path,ssd_rate,loss_rate)

````


**root(species_tree,gene_tree)** | This function accepts a gene tree and species tree. This function will pick the best rooting of the gene tree to minimize the number of loss and duplication events.

```python
species_tree = PhyloTree("example_species_tree.newick")

gene_tree = PhyloTree("example_gene_tree.newick")

rooted_gene_tree = wgdtree.root(species_tree,gene_tree)
````

**place_wgd(species_tree, gene_tree)** |  This function accepts a rooted gene tree and a species tree labelled with one pair of consecutive whole genome duplication events. The function will return a copy of the gene tree with labels added to nodes that represent wgd events. The added events will either be marked with a 'P' for present or 'M' for missing. All annotations will be made in NHX format.

```python
species_tree = PhyloTree("example_species_tree.newick")

gene_tree = PhyloTree("example_gene_tree.newick")

labelled_gene_tree = wgdtree.place_wgd(species_tree, gene_tree)
````


**rrates(labelled_gene_tree, pSpecies)** |  This function accepts a gene tree labelled using the place_wgd() method and a string corresponding to a taxa on the given gene tree. This function then returns two floating point values representing the conditional probably of genes being retained retained and lossed retained of the given taxa. 


```python

pSpecies = "Dinosourous_rex"

results = wgdtree.rrates(labeled_gene_tree, pSpecies)

retained_retained = results[0]
                
lossed_retained = results[1]
````

***run(list_of_gene_tree,species_tree)*** | This function accepts an array of gene tree files in newick format and a species tree 
labelled with all whole genome duplication events on the lineage. The function will calculate the conditional probability statistic 
for all consecutive pairs of events on the species tree. The function will return a dictionary where species names as they appear on the
species trees are used as keys that correspond to a list of tuples of the form (P(L|R),P(R|R)). There will be a seperate tuple for each event 
pair that species was present for.

```python

list_of_gene_trees = ["gene_tree0.newick","gene_tree1.newick","gene_tree2.newick","gene_tree3.newick"]

species_tree = PhyloTree("species_tree")

results_dict = wgdtree.run(list_of_gene_trees,species_tree)
```

example results_dict.txt
```
{

'species_A': [], 

'species_B': [(0.647454313863336, 0.5535554282888006)], 

'species_C': [(0.448686400377538, 0.2861785837493464)], 

'species_D': [(0.37466663338849644, 0.3420519777168986)], 

'species_E': [(0.205555576769, 0.2314676413445378)], 

'species_F': [(0.21592868256170894, 0.18623158187378372), (0.3076854227579211, 0.2671778530173272), (0.29352608422375864, 0.2732342007434944)], 

'species_G': [(0.26377679954349065, 0.23002217890987298), (0.3769475818129163, 0.33331715367440057), (0.36203343256239984, 0.33394428152492667)]

}

```


*all trees must be in .newick format, and the ete3 package is required.

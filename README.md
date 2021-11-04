# wgdtree Package


install with ```$ pip install wgdtree ```

to use the package in your code add import wgdtree

This package provides the following functions:


**find_wgd(species_tree, gene_tree)** |  This function accepts a species tree and a rooted gene tree. The function will return a copy of the gene tree 
with labels added to nodes that represent wgd event. The added events will either be marked with a 'P' for present or 'M' for missing. It is important that the species tree have labels on branches were the wgd events occured.

```python
find_wgd(species_tree, gene_tree)
````



**sim(species_tree,n,path,ssd_rate,loss_rate)** | This function accepts a species tree, an integer value 'n', a path, a loss rate floating point and a ssd rate also floating point. The function will simulate n statistically probable gene tree and then write them as newick files to the given path. This function allows for wgd duplications that are labelled on the species tree to be accounted for.


```python
sim(species_tree,n,path,ssd_rate,loss_rate)
````



**root(species_tree,gene_tree,mode)** | This function accepts a gene tree and species tree. This function will pick the best rooting of the gene tree to minimize the number of loss and duplication events. Sudo code for this function was taken from [xx].  


```python
root(species_tree,gene_tree,mode)
````


**rrates(labeled_gene_tree)** |  This function accepts a gene tree labelled using the find_wgd() method. This function then returns two floating point values representing the conditional probably of retained retained and lossed retained. 


```python
rrates(labeled_gene_tree)
````



*all trees must be in .newick format, and the ete3 package is required.

import sys
import os
from species import species
from Bio import Phylo
M = []
O = []
infile = sys.argv[1]
tree = Phylo.read(infile, "newick")
a = tree.get_terminals()
for g in a:
    G = int(g.name)
    gname = species[G + 1].strip()
    ggname = gname.replace(" ", "")
    Gn = ggname + "_" + g.name
    g.name = Gn
    Phylo.write([tree], "labeled" + infile, "newick")

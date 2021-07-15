from wgdtree import root, rrates, place_wgd
from ete3 import PhyloTree, Tree
import sys

g_file = sys.argv[1]
s_file = sys.argv[2]


g_tree = PhyloTree(g_file, format=1)
s_tree = PhyloTree(s_file, format=1)

g_tree = root(g_tree,s_tree)

place_wgd(s_tree,g_tree)

print(rrates(g_tree,"P.halli"))

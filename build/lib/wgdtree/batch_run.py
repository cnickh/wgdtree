from ete3 import Tree
from wgdtree.retention_rates import rrates, place_wgd

def run(list_of_gene_trees, species_tree):

    results_dic = {}
    paired_trees = []
    pSpecies = []

    for l in species_tree.get_leaves():
        results_dic[l.name] = []

    clean_tree = species_tree.copy()

    for n in clean_tree.iter_descendants():
        if('WGD' in n.features):
            n.del_feature('WGD')

    i = 0
    for n in species_tree.traverse():

        if('WGD' in n.features and n.WGD == 'Y'):
            for x in n.iter_descendants():

                if('WGD' in x.features and x.WGD == 'Y'):

                    tmp_tree = clean_tree.copy()

                    n1 = tmp_tree&n.name

                    n2 = tmp_tree&x.name

                    n1.add_feature('WGD', 'Y')

                    n2.add_feature('WGD', 'Y')

                    paired_trees.append(tmp_tree)

                    i+=1
                    leaves = []

                    for l in x.get_leaves():
                        leaves.append(l.name.split('_')[0])

                    pSpecies.append(leaves)

    i=0
    for t in paired_trees:  #for each pair of events
        for s in pSpecies[i]: #for each species present for the pair
            l_poss = 0
            r_poss = 0
            l = 0
            r = 0
            q = 0

            for g in list_of_gene_trees: #for each gene tree
                tmp_gene = g.copy()
                place_wgd(t,tmp_gene)

                results = rrates(tmp_gene,s)

                q+=1
                if(results != -1):
                    l_poss += results[1][1]
                    r_poss += results[0][1]
                    l += results[1][0]
                    r += results[0][0]

            if(l_poss == 0 and r_poss == 0):
                results_dic[s].append((0,0))
            elif (l_poss == 0):
                results_dic[s].append((0,r/r_poss))
            elif (r_poss == 0):
                results_dic[s].append((l/l_poss,0))
            else:
                results_dic[s].append((l/l_poss,r/r_poss))



        i+=1

    return results_dic

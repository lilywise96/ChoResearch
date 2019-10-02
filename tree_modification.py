"""
Name: tree_modification.py
Author: Lily Wise


"""


# Swaps a dictionary that has a key and a value where the value is a list
# or a set.
#
# param: key_to_value - a dictionary; key: key, value: set or list of values
# return: value_to_key - a dictionary; key: previous value, value: set or list of previous keys
def swap_key_value(key_to_value):
    value_to_key = {}

    for key in key_to_value:
        for value in key_to_value[key]:
            if value not in value_to_key:
                value_to_key[value] = set()
            value_to_key[value].add(key)

    return value_to_key


# Find all leaf nodes and add all of the genes to their parents.
def gene_to_all_parents(tree_child_parent, gene_to_terms):
    tree_parent_child = swap_key_value(tree_child_parent)
    all_nodes = set()
    terms_to_genes = swap_key_value(gene_to_terms)
    leaf_nodes = set()
    to_check = set()

    for node in set(tree_parent_child.keys()):
        all_nodes.add(node)
    for node in set(tree_child_parent.keys()):
        all_nodes.add(node)

    for node in all_nodes:
        if node not in set(tree_parent_child.keys()):
            leaf_nodes.add(node)

    print(leaf_nodes)

    for node in leaf_nodes:
        for parent in tree_child_parent[node]:
            to_check.add(parent)

    while len(to_check) != 0:
        checking = to_check.pop()
        print(checking)
        if checking in tree_child_parent:
            for parent in tree_child_parent[checking]:
                to_check.add(parent)
            terms_to_genes = gene_to_parent(checking, tree_parent_child, terms_to_genes)

    return terms_to_genes


# Add all gene_carry genes to the current node. Add current nodes
# genes to gene_carry and then call for all parents of current node
# this function until there are no parents left to call.
#
# param: node - the current node in the tree
# param: tree - a dictionary; key: term, value: set of parent terms
# param: gene_to_terms - a dictionary; key: gene, value: set of terms that it is annotated to
# return: new_gene_to_terms - updated gene_to_terms with the node that was called having all of its children's genes
def gene_to_parent(node, tree, terms_to_genes):
    new_terms_to_genes = terms_to_genes
    print(new_terms_to_genes)

    for child in tree[node]:
        if node not in new_terms_to_genes:
            new_terms_to_genes[node] = set()
        if child in new_terms_to_genes:
            for c in new_terms_to_genes[child]:
                new_terms_to_genes[node].add(c)

    return new_terms_to_genes

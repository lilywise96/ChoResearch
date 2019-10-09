"""
Name: tree_modification.py
Author: Lily Wise

All functions that modify or determine information about the annotation or ontology trees.
"""

from math import log10


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


# Joins all trees and switches from term to genes into gene to terms.
#
# param: hp_tg - the human phenotype tree, term to genes
# param: bp_tg - the biological process tree, term to genes
# param: mf_tg - the molecular function tree, term to genes
#
# returns: the joined tree, gene to terms
def join_gt(hp_tg, bp_tg, mf_tg):
    hp_gt = swap_key_value(hp_tg)
    bp_gt = swap_key_value(bp_tg)
    mf_gt = swap_key_value(mf_tg)

    all_gt = {}
    for gene in hp_gt:
        if gene not in all_gt:
            all_gt[gene] = set()
        for term in hp_gt[gene]:
            all_gt[gene].add(term)

    for gene in bp_gt:
        if gene not in all_gt:
            all_gt[gene] = set()
        for term in bp_gt[gene]:
            all_gt[gene].add(term)

    for gene in mf_gt:
        if gene not in all_gt:
            all_gt[gene] = set()
        for term in mf_gt[gene]:
            all_gt[gene].add(term)

    return all_gt


# Find all leaf nodes and add all of the genes to their parents.
#
# param: tree_child_parent - tree that is child to parent
# param: gene_to_terms - gene to terms it is annotated to
# returns: terms to genes tree with all ancestors having all the terms below theirs genes
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

    for node in leaf_nodes:
        for parent in tree_child_parent[node]:
            to_check.add(parent)

    while len(to_check) != 0:
        checking = to_check.pop()
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

    for child in tree[node]:
        if node not in new_terms_to_genes:
            new_terms_to_genes[node] = set()
        if child in new_terms_to_genes:
            for c in new_terms_to_genes[child]:
                new_terms_to_genes[node].add(c)

    return new_terms_to_genes


# Calculates information content of every term of the tree
#
# param: tree - the tree; key: term, value: genes annotated to that term (with transitive property)
# return: term_ic - dictionary; key: term, value: information content of that term
def calculate_ic(tree):
    term_ic = {}  # key: term, value: ic
    distinct_genes = set()

    for term in tree:
        term_ic[term] = len(tree[term])
        distinct_genes.add(tree[term])

    for term in term_ic:
        term_ic[term] /= len(distinct_genes)
        term_ic[term] = -log10(term_ic[term])

    return term_ic

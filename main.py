"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases.
"""


from ontology_parsing import hpo_parsing_onto, parsing_go, testing_ontology_parsing
from annotation_parsing import hpo_parsing_ann, parsing_ann, testing_annotation_parsing
from tree_modification import swap_key_value, gene_to_all_parents, join_gt
from apriori_algorithm import apriori
from association_creation import create_associations
from math import ceil
import sys


def create_onto_ann(gene_term_output_filename):
    # Read in ontologies. terms to parents
    hpo_terms_parents = hpo_parsing_onto("./input/hp.obo.txt")
    bp_terms_parents, mf_terms_parents, cc_terms_parents = parsing_go("./input/go.obo")
    # testing_terms_parents = testing_ontology_parsing("./input/testing_ontology")

    # Read in annotations.
    hp_gt = hpo_parsing_ann("./input/hpo_genes_to_phenotype.txt")
    gene_syn, bp_gt, mf_gt, cc_gt = parsing_ann("./input/goa_human.gaf")
    # testing_gt = testing_annotation_parsing("./input/testing_annotation")

    # Recursively add genes to parents.
    hp_tg = gene_to_all_parents(hpo_terms_parents, hp_gt)
    bp_tg = gene_to_all_parents(bp_terms_parents, bp_gt)
    mf_tg = gene_to_all_parents(mf_terms_parents, mf_gt)
    # cc_tg = gene_to_all_parents(cc_terms_parents, cc_gt)
    # testing_tg = gene_to_all_parents(testing_terms_parents, testing_gt)

    # Save terms for each.
    hp_terms = set(hp_tg.keys())
    bp_terms = set(bp_tg.keys())
    mf_terms = set(mf_tg.keys())

    # Join all term to gene and make them gene to term.
    all_gt = join_gt(hp_tg, bp_tg, mf_tg)

    output_file = open(gene_term_output_filename, "w")
    for gene in all_gt:
        output_file.write(gene)
        for term in all_gt[gene]:
            output_file.write("\t")
            output_file.write(term)
        output_file.write("\n")
    output_file.close()
    return all_gt


def read_onto_ann(filename):
    file = open(filename, "r")
    gt = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        gt[cols[0]] = set()
        for i in range(1, len(cols)):
            gt[cols[0]].add(cols[i])
    file.close()
    return gt


def create_freq_itemsets(filename, all_gt):
    all_terms = set()
    for gene in all_gt:
        for term in all_gt[gene]:
            all_terms.add(term)

    freq_itemsets = apriori(all_gt, all_terms, min_support)

    file = open("freq_itemsets.txt", "w")
    for itemset in freq_itemsets:
        for item in range(0, len(itemset)):
            if item != 0:
                file.write("\t")
            file.write(itemset[item])
        file.write("\n")
    file.close()

    return freq_itemsets


def read_freq_itemsets(filename):
    freq_itemsets = []
    file = open(filename, "r")
    for line in file:
        items = line.split("\t")
        items[len(items) - 1] = items[len(items) - 1][0:-1]

        itemset = set()
        for item in items:
            itemset.add(item)
        freq_itemsets.append(itemset)
    file.close()

    return freq_itemsets


def create_new_associations(all_gt, freq_itemsets, min_confidence, filename):
    final_associations = create_associations(all_gt, freq_itemsets, min_confidence)
    file = open(filename, "w")

    for association in final_associations:
        for associate in range(0, len(association)):
            if associate != 0:
                file.write("\t")
            file.write(association[associate])
        file.write("\n")
    file.close()

    return final_associations


def read_associations(filename):
    final_associations = []
    file = open(filename, "r")
    for line in file:
        association = line.split("\t")
        association[len(association) - 1] = association[len(association) - 1][0:-1]

        final_association = set()
        for assoc in association:
            final_association.add(assoc)
        final_associations.append(final_association)
    file.close()

    return final_associations


if len(sys.argv) == 7:
    print("Correct number of variables.")

    recreate_onto_ann = sys.argv[1]
    recreate_freq_itemsets = sys.argv[2]
    recreate_associations = sys.argv[3]
    min_support = float(sys.argv[4])
    min_confidence = float(sys.argv[5])
    min_information_content = float(sys.argv[6])

    if recreate_onto_ann == "true":
        all_gt = create_onto_ann("gene_term.txt")
    else:
        all_gt = read_onto_ann("gene_term.txt")

    if recreate_freq_itemsets == "true":
        freq_itemsets = create_freq_itemsets("freq_itemsets.txt", all_gt)
    else:
        freq_itemsets = read_freq_itemsets("freq_itemsets.txt")

    if recreate_associations == "true":
        final_associations = create_new_associations(all_gt, freq_itemsets, min_confidence, "associations.txt")
    else:
        final_associations = read_associations("associations.txt")

    print("Done")

else:
    print("Not the correct number of variables: (rewrite_onto_ann, rewrite_freq_itemsets, "
          "min_support, min_confidence, min_information_content)")

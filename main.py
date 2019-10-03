"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases.
"""


from ontology_parsing import hpo_parsing_onto, parsing_go, testing_ontology_parsing
from annotation_parsing import hpo_parsing_ann, parsing_ann, testing_annotation_parsing
from tree_modification import swap_key_value, gene_to_all_parents, join_gt
from apriori_algorithm import apriori
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

    return gt


if len(sys.argv) == 5:
    print("Correct number of variables.")

    recreate = sys.argv[1]
    min_support = float(sys.argv[2])
    min_confidence = float(sys.argv[3])
    min_information_content = float(sys.argv[4])

    if recreate == "true":
        all_gt = create_onto_ann("gene_term.txt")
    else:
        all_gt = read_onto_ann("gene_term.txt")

    print("Calculating frequent itemsets.")
    freq_itemsets = apriori(all_gt, set(all_gt.keys()), min_support)
    # print(freq_itemsets)

else:
    print("Not the correct number of variables: (rewrite, min_support, min_confidence, min_information_content)")

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

if len(sys.argv) == 4:
    print("Correct number of variables.")

    min_support = float(sys.argv[1])
    min_confidence = float(sys.argv[2])
    min_information_content = float(sys.argv[3])

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

    print("Calculating frequent itemsets.")
    # print(set(all_gt.keys()))
    # print(min_support)
    freq_itemsets = apriori(all_gt, set(all_gt.keys()), min_support)
    # print(freq_itemsets)

else:
    print("Not the correct number of variables: (min_support, min_confidence, min_information_content)")

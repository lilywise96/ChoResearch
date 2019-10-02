"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases.
"""


from ontology_parsing import hpo_parsing_onto, parsing_go, testing_ontology_parsing
from annotation_parsing import hpo_parsing_ann, parsing_ann, testing_annotation_parsing
from tree_modification import swap_key_value, gene_to_all_parents
from math import ceil
import sys

if len(sys.argv) == 4:
    print("Correct number of variables.")

    min_support = sys.argv[1]
    min_confidence = sys.argv[2]
    min_information_content = sys.argv[3]

    # Read in ontologies. terms to parents
    # hpo_terms_parents = hpo_parsing_onto("hp.obo.txt")
    # go_terms_parents = parsing_go("go.obo")
    testing_terms_parents = testing_ontology_parsing("testing_ontology")

    # Read in annotations.
    # hp_gt = hpo_parsing_ann("hpo_genes_to_phenotype.txt")
    # gene_syn, bp_gt, mf_gt, cc_gt = parsing_ann("goa_human.gaf")
    testing_gt = testing_annotation_parsing("testing_annotation")

    # Recursively add genes to parents.
    testing_tg = gene_to_all_parents(testing_terms_parents, testing_gt)


else:
    print("Not the correct number of variables: (min_support, min_confidence, min_information_content)")

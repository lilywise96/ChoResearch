"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases.
"""


from ontology_parsing import hpo_parsing_onto, parsing_go
from annotation_parsing import hpo_parsing_ann, parsing_ann
from math import ceil
import sys

if len(sys.argv) == 4:
    print("Correct number of variables.")

    min_support = sys.argv[1]
    min_confidence = sys.argv[2]
    min_information_content = sys.argv[3]

    # Read in ontologies.
    hpo_terms_parents = hpo_parsing_onto("hp.obo.txt")
    go_terms_parents = parsing_go("go.obo")

    # Read in annotations.
    hp_gt = hpo_parsing_ann("hpo_genes_to_phenotype.txt")
    gene_syn, bp_gt, mf_gt, cc_gt = parsing_ann("goa_human.gaf")

else:
    print("Not the correct number of variables: (min_support, min_confidence, min_information_content)")

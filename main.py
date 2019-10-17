"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases. Should calculate associations with
support 4% - 10%, coverage 4% - 10%, and confidence 20% - 50%.
"""

# BP -> BP and MF -> MF and HPO -> HPO and BP -> HPO and MF -> HPO

from ontology_parsing import hpo_parsing_onto, parsing_go, testing_ontology_parsing
from annotation_parsing import hpo_parsing_ann, parsing_ann, testing_annotation_parsing
from tree_modification import gene_to_all_parents, join_gt, calculate_ic, terms_to_all_parents, all_specificity
from apriori_algorithm import apriori
from association_creation import create_associations
from math import ceil
import sys

# Directories
input_direct = "./input/"
created_direct = "./created/"

# File names
# Specificity Storage
hp_spec_filename = created_direct + "hp_spec.txt"
bp_spec_filename = created_direct + "bp_spec.txt"
mf_spec_filename = created_direct + "mf_spec.txt"

# Information Content Storage
hp_ic_filename = created_direct + "hp_ic.txt"
bp_ic_filename = created_direct + "bp_ic.txt"
mf_ic_filename = created_direct + "mf_ic.txt"

# Transitivity of Terms Storage
hp_trans_filename = created_direct + "hp_trans.txt"
bp_trans_filename = created_direct + "bp_trans.txt"
mf_trans_filename = created_direct + "mf_trans.txt"

# Gene to Terms All
gene_term_filename = created_direct + "gene_term.txt"

# Ontology Given
hp_ontology_filename = input_direct + "hp.obo.txt"
g_ontology_filename = input_direct + "go.obo"

# Annotation Given
hp_annotations_filename = input_direct + "hpo_genes_to_phenotype.txt"
g_annotations_filename = input_direct + "goa_human.gaf"

# Modified with Program Parameters
# Frequent Itemsets
freq_itemsets_filename = created_direct + "freq_itemsets_1.txt"

# Associations
associations_filename = created_direct + "associations_1.txt"


# Creates ontology and writes it to an output file, as well as calculates information content.
#
# param: gene_term_output_filename - the file the ontology is written to
# return: all_gt - dictionary; key: gene, value: set of terms
def create_onto_ann():
    # Read in ontologies. terms to parents
    hpo_terms_parents = hpo_parsing_onto(hp_ontology_filename)
    bp_terms_parents, mf_terms_parents, cc_terms_parents = parsing_go(g_ontology_filename)

    # Read in annotations.
    hp_gt = hpo_parsing_ann(hp_annotations_filename)
    gene_syn, bp_gt, mf_gt, cc_gt = parsing_ann(g_annotations_filename)

    # Creates transitive trees
    hpo_terms_parents_trans = terms_to_all_parents(hpo_terms_parents)
    bp_terms_parents_trans = terms_to_all_parents(bp_terms_parents)
    mf_terms_parents_trans = terms_to_all_parents(mf_terms_parents)

    # Recursively add genes to parents.
    hp_tg = gene_to_all_parents(hpo_terms_parents, hp_gt)
    bp_tg = gene_to_all_parents(bp_terms_parents, bp_gt)
    mf_tg = gene_to_all_parents(mf_terms_parents, mf_gt)

    # Calculate Information Content
    hp_ic = calculate_ic(hp_tg)
    bp_ic = calculate_ic(bp_tg)
    mf_ic = calculate_ic(mf_tg)

    # Save terms for each
    hp_terms = set(hp_tg.keys())
    bp_terms = set(bp_tg.keys())
    mf_terms = set(mf_tg.keys())

    # Specificity of all nodes
    hp_spec = all_specificity(hpo_terms_parents_trans, hp_ic)
    bp_spec = all_specificity(bp_terms_parents_trans, bp_ic)
    mf_spec = all_specificity(mf_terms_parents_trans, mf_ic)

    # Join all term to gene and make them gene to term.
    all_gt = join_gt(hp_tg, bp_tg, mf_tg)

    output_file = open(gene_term_filename, "w")
    for gene in all_gt:
        output_file.write(gene)
        for term in all_gt[gene]:
            output_file.write("\t")
            output_file.write(term)
        output_file.write("\n")
    output_file.close()

    # Write specificity to files.
    output_file = open(hp_spec_filename, "w")
    for term in hp_spec:
        output_file.write(term+"\t"+str(hp_spec[term])+"\n")
    output_file.close()

    output_file = open(bp_spec_filename, "w")
    for term in bp_spec:
        output_file.write(term+"\t"+str(bp_spec[term])+"\n")
    output_file.close()

    output_file = open(mf_spec_filename, "w")
    for term in mf_spec:
        output_file.write(term+"\t"+str(mf_spec[term])+"\n")
    output_file.close()

    # Write information content to files.
    output_file = open(hp_ic_filename, "w")
    for term in hp_ic:
        output_file.write(term + "\t" + str(hp_ic[term]) + "\n")
    output_file.close()

    output_file = open(bp_ic_filename, "w")
    for term in bp_ic:
        output_file.write(term + "\t" + str(bp_ic[term]) + "\n")
    output_file.close()

    output_file = open(mf_ic_filename, "w")
    for term in mf_ic:
        output_file.write(term + "\t" + str(mf_ic[term]) + "\n")
    output_file.close()

    all_ic = {}
    for term in hp_ic:
        all_ic[term] = hp_ic[term]
    for term in bp_ic:
        all_ic[term] = bp_ic[term]
    for term in mf_ic:
        all_ic[term] = mf_ic[term]

    all_spec = {}
    for term in hp_spec:
        all_spec[term] = hp_spec[term]
    for term in bp_spec:
        all_spec[term] = bp_spec[term]
    for term in mf_spec:
        all_spec[term] = mf_spec[term]

    return all_gt, all_spec, all_ic


def read_onto_ann():
    file = open(gene_term_filename, "r")
    gt = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        gt[cols[0]] = set()
        for i in range(1, len(cols)):
            gt[cols[0]].add(cols[i])
    file.close()

    # Read Specificity Files
    file = open(hp_spec_filename, "r")
    hp_spec = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        hp_spec[cols[0]] = float(cols[1])
    file.close()

    file = open(bp_spec_filename, "r")
    bp_spec = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        bp_spec[cols[0]] = float(cols[1])
    file.close()

    file = open(mf_spec_filename, "r")
    mf_spec = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        mf_spec[cols[0]] = float(cols[1])
    file.close()

    # Read Information Content Files
    file = open(hp_ic_filename, "r")
    hp_ic = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        hp_ic[cols[0]] = float(cols[1])
    file.close()

    file = open(bp_ic_filename, "r")
    bp_ic = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        bp_ic[cols[0]] = float(cols[1])
    file.close()

    file = open(mf_ic_filename, "r")
    mf_ic = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        mf_ic[cols[0]] = float(cols[1])
    file.close()

    all_ic = {}
    for term in hp_ic:
        all_ic[term] = hp_ic[term]
    for term in bp_ic:
        all_ic[term] = bp_ic[term]
    for term in mf_ic:
        all_ic[term] = mf_ic[term]

    all_spec = {}
    for term in hp_spec:
        all_spec[term] = hp_spec[term]
    for term in bp_spec:
        all_spec[term] = bp_spec[term]
    for term in mf_spec:
        all_spec[term] = mf_spec[term]

    return all_gt, all_spec, all_ic


# Create frequent itemsets.
def create_freq_itemsets(filename, all_gt, min_support, min_weighted_support,
                         min_information_content, all_spec, all_ic):
    all_terms = set()
    for gene in all_gt:
        for term in all_gt[gene]:
            all_terms.add(term)

    freq_itemsets = apriori(all_gt, all_terms, min_support, min_weighted_support,
                            min_information_content, all_spec, all_ic)

    file = open(filename, "w")
    file.write("Min Support - "+str(min_support)+"\n")
    file.write("Min Information Content - " + str(min_information_content) + "\n")
    file.write("Min Weighted Support - " + str(min_weighted_support) + "\n")
    for itemset in freq_itemsets:
        for item in range(0, len(itemset)):
            if item != 0:
                file.write("\t")
            file.write(itemset[item])
        file.write("\n")
    file.close()

    return freq_itemsets


# Read frequent itemsets.
def read_freq_itemsets(filename):
    freq_itemsets = []
    file = open(filename, "r")
    count = 3
    for line in file:
        if count > 0:
            count -= 1
        else:
            items = line.split("\t")
            items[len(items) - 1] = items[len(items) - 1][0:-1]

            itemset = set()
            for item in items:
                itemset.add(item)
            freq_itemsets.append(itemset)
    file.close()

    return freq_itemsets


def create_new_associations(all_gt, freq_itemsets, min_confidence, min_coverage, filename, all_spec):
    final_associations = create_associations(all_gt, freq_itemsets, min_confidence, min_coverage, all_spec)
    file = open(filename, "w")

    file.write("Min Coverage - " + str(min_coverage) + "\n")
    file.write("Min Confidence - " + str(min_confidence) + "\n")

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

    count = 2
    for line in file:
        if count > 0:
            count -= 1
        else:
            association = line.split("\t")
            association[len(association) - 1] = association[len(association) - 1][0:-1]

            final_association = set()
            for assoc in association:
                final_association.add(assoc)
            final_associations.append(final_association)
    file.close()

    return final_associations


if len(sys.argv) == 9:
    print("Correct number of variables.")

    recreate_onto_ann = sys.argv[1]
    recreate_freq_itemsets = sys.argv[2]
    recreate_associations = sys.argv[3]
    min_support = float(sys.argv[4])
    min_weighted_support = float(sys.argv[5])
    min_confidence = float(sys.argv[6])
    min_information_content = float(sys.argv[7])
    min_coverage = float(sys.argv[8])

    print("Min_Support: "+str(min_support))
    print("Min_Weighted_Support: "+str(min_weighted_support))
    print("Min_Confidence: "+str(min_confidence))
    print("Min_Info_Content: "+str(min_information_content))
    print("Min_Coverage: "+str(min_coverage))

    # freq_itemsets_filename += str(ceil(min_support*100)) + "_" + str(ceil(min_weighted_support*100)) + "_" + \
    #                           str(ceil(min_information_content*100)) + ".txt"
    # associations_filename += str(ceil(min_confidence*100)) + "_" + str(ceil(min_coverage*100)) + ".txt"

    if recreate_onto_ann == "true":
        all_gt, all_spec, all_ic = create_onto_ann()
    else:
        all_gt, all_spec, all_ic = read_onto_ann()

    if recreate_freq_itemsets == "true":
        freq_itemsets = create_freq_itemsets(freq_itemsets_filename, all_gt, min_support, min_weighted_support,
                                             min_information_content, all_spec, all_ic)
    else:
        freq_itemsets = read_freq_itemsets(freq_itemsets_filename)

    min_coverage_count = ceil(min_coverage * len(freq_itemsets))
    min_confidence_count = ceil(min_confidence * len(freq_itemsets))

    if recreate_associations == "true":
        final_associations = create_new_associations(all_gt, freq_itemsets, min_confidence, min_coverage,
                                                     associations_filename, all_spec)
    else:
        final_associations = read_associations(associations_filename)

    print("Done")

else:
    print("Not the correct number of variables: (rewrite_onto_ann, rewrite_freq_itemsets, rewrite_associations, "
          "min_support, min_weighted_support, min_confidence, min_information_content, min_coverage)")

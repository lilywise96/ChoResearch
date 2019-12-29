"""
Filename: main.py
Author: Lily Wise

Calls other functions to find associations between genes and terms for diseases. Should calculate associations with
support 4% - 10%, coverage 4% - 10%, and confidence 20% - 50%.
"""

# BP -> BP and MF -> MF and HPO -> HPO and BP -> HPO and MF -> HPO
# Accuracy Measuring?

from ontology_parsing import hpo_parsing_onto, parsing_go, testing_ontology_parsing
from annotation_parsing import hpo_parsing_ann, parsing_ann, testing_annotation_parsing
from tree_modification import gene_to_all_parents, join_gt, calculate_ic, terms_to_all_parents, \
    all_specificity, swap_key_value, join_two
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
gene_term_bp_filename = created_direct + "gene_term_bp.txt"
gene_term_mf_filename = created_direct + "gene_term_mf.txt"
gene_term_hp_filename = created_direct + "gene_term_hp.txt"

# Ontology Given
hp_ontology_filename = input_direct + "hp.obo.txt"
g_ontology_filename = input_direct + "go.obo"

# Annotation Given
hp_annotations_filename = input_direct + "hpo_genes_to_phenotype.txt"
g_annotations_filename = input_direct + "goa_human.gaf"


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
    hp_gt = swap_key_value(hp_tg)
    bp_gt = swap_key_value(bp_tg)
    mf_gt = swap_key_value(mf_tg)

    output_file = open(gene_term_filename, "w")
    for gene in all_gt:
        output_file.write(gene)
        for term in all_gt[gene]:
            output_file.write("\t")
            output_file.write(term)
        output_file.write("\n")
    output_file.close()

    output_file = open(gene_term_bp_filename, "w")
    for gene in bp_gt:
        output_file.write(gene)
        for term in bp_gt[gene]:
            output_file.write("\t")
            output_file.write(term)
        output_file.write("\n")
    output_file.close()

    output_file = open(gene_term_mf_filename, "w")
    for gene in mf_gt:
        output_file.write(gene)
        for term in mf_gt[gene]:
            output_file.write("\t")
            output_file.write(term)
        output_file.write("\n")
    output_file.close()

    output_file = open(gene_term_hp_filename, "w")
    for gene in hp_gt:
        output_file.write(gene)
        for term in hp_gt[gene]:
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

    return all_gt, all_spec, all_ic, bp_gt, mf_gt, hp_gt


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

    file = open(gene_term_bp_filename, "r")
    bp_gt = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        bp_gt[cols[0]] = set()
        for i in range(1, len(cols)):
            bp_gt[cols[0]].add(cols[i])
    file.close()

    file = open(gene_term_mf_filename, "r")
    mf_gt = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        mf_gt[cols[0]] = set()
        for i in range(1, len(cols)):
            mf_gt[cols[0]].add(cols[i])
    file.close()

    file = open(gene_term_hp_filename, "r")
    hp_gt = {}
    for line in file:
        cols = line.split("\t")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]

        hp_gt[cols[0]] = set()
        for i in range(1, len(cols)):
            hp_gt[cols[0]].add(cols[i])
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

    return gt, all_spec, all_ic, bp_gt, mf_gt, hp_gt


# Create frequent itemsets.
def create_freq_itemsets(filename, possible_left, all_gt, min_support, min_weighted_support,
                         min_information_content, all_spec, all_ic):
    all_terms = set()
    for gene in all_gt:
        for term in all_gt[gene]:
            all_terms.add(term)

    print("ALL GT: ", end="")
    print(all_gt)

    freq_itemsets = apriori(all_gt, all_terms, min_support, min_weighted_support,
                            min_information_content, all_spec, all_ic)

    all_itemsets = []
    for size_freq in freq_itemsets:
        for itemset in freq_itemsets[size_freq]:
            all_itemsets.append(itemset)

    freq_itemsets = all_itemsets

    for itemset in freq_itemsets:
        found = False
        for item in itemset:
            if item in possible_left:
                found = True

        if not found:
            freq_itemsets.remove(itemset)

    file = open(filename, "w")
    file.write("Min Support - "+str(min_support)+"\n")
    file.write("Min Information Content - " + str(min_information_content) + "\n")
    file.write("Min Weighted Support - " + str(min_weighted_support) + "\n")
    for itemset in freq_itemsets:
        item_list = list(itemset)
        for i in range(0, len(item_list)):
            if i != 0:
                file.write("\t")
            file.write(item_list[i])
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


def create_new_associations(left_terms, right_terms, all_gt, freq_itemsets, min_confidence, min_coverage,
                            filename, all_spec):
    final_associations = create_associations(left_terms, right_terms, all_gt, freq_itemsets, min_confidence,
                                             min_coverage, all_spec)
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


def general_main(freq_file_ext, association_file_ext, recreate_onto_ann, recreate_freq_itemsets, tree,
                 min_support, min_weighted_support, min_confidence, min_information_content, min_coverage):

    freq_itemsets_filename = created_direct + "freq_itemsets_" + str(freq_file_ext) + ".txt"
    associations_filename = created_direct + "associations_" + str(association_file_ext) + ".txt"
    information_filename = "info_on_files.txt"

    info_file = open(information_filename, "a+")

    print("Frequent itemsets filename: "+str(freq_itemsets_filename))
    print("Associations filename: "+str(associations_filename))
    print("Tree: "+str(tree))
    print("Minimum support: "+str(min_support))
    print("Minimum weighted support: "+str(min_weighted_support))
    print("Minimum confidence: "+str(min_confidence))
    print("Minimum information content: "+str(min_information_content))
    print("Minimum coverage: "+str(min_coverage))

    info_file.write("Frequent itemsets filename: " + str(freq_itemsets_filename) + "\n")
    info_file.write("Associations filename: " + str(associations_filename) + "\n")
    info_file.write("Tree: " + str(tree) + "\n")
    info_file.write("Minimum support: " + str(min_support) + "\n")
    info_file.write("Minimum weighted support: " + str(min_weighted_support) + "\n")
    info_file.write("Minimum confidence: " + str(min_confidence) + "\n")
    info_file.write("Minimum information content: " + str(min_information_content) + "\n")
    info_file.write("Minimum coverage: " + str(min_coverage) + "\n")
    info_file.write("\n\n")

    info_file.close()

    if recreate_onto_ann == "true":
        all_gt, all_spec, all_ic, bp_gt, mf_gt, hp_gt = create_onto_ann()
    else:
        all_gt, all_spec, all_ic, bp_gt, mf_gt, hp_gt = read_onto_ann()

    possible_left = set()
    possible_right = set()
    if tree == 'bp':
        for gene in bp_gt:
            possible_left = bp_gt[gene].union(possible_left)
        for gene in hp_gt:
            possible_right = hp_gt[gene].union(possible_right)
    elif tree == 'mf':
        for gene in mf_gt:
            possible_left = mf_gt[gene].union(possible_left)
        for gene in hp_gt:
            possible_right = hp_gt[gene].union(possible_right)
    elif tree == 'hp':
        for gene in hp_gt:
            possible_left = hp_gt[gene].union(possible_left)
        possible_right = possible_left
    else:
        for gene in all_gt:
            possible_left = all_gt[gene].union(possible_left)
        possible_right = possible_left

    if recreate_freq_itemsets == "true":
        if tree == 'bp':
            all_gt = join_two(bp_gt, hp_gt)
        elif tree == 'mf':
            all_gt = join_two(mf_gt, hp_gt)
        elif tree == 'hp':
            all_gt = hp_gt

        freq_itemsets = create_freq_itemsets(freq_itemsets_filename, possible_left, all_gt,
                                             min_support, min_weighted_support, min_information_content,
                                             all_spec, all_ic)
    else:
        freq_itemsets = read_freq_itemsets(freq_itemsets_filename)

    create_new_associations(possible_left, possible_right, all_gt, freq_itemsets, min_confidence,
                            min_coverage, associations_filename, all_spec)

    print("Done")

"""
File: annotation_parsing.py
Author: Lily Wise, Joseph Chang

This file parses annotation files for hpo and for go.
"""


# This function parses the hpo annotation file. It pulls the gene annotation
# and terms it is associated to.
#
# param: filename - the hpo annotation file
# return: gene_term_id - a dictionary; key: gene, value: array of terms
def hpo_parsing_ann(filename):
    file = open(filename, "r")
    gene_id_symbol = {}
    gene_term_id = {}

    for line in file:
        if not line.startswith('#'):
            columns = line.strip().split('\t')
            gene_id = columns[0]
            gene_symbol = columns[1]
            term_id = columns[3][0:10]

            if gene_id not in gene_id_symbol.keys():
                gene_id_symbol[gene_id] = gene_symbol
                gene_term_id[gene_symbol] = []

            gene_term_id[gene_symbol].append(term_id)

    return gene_term_id


# This function parses the gene annotation file. It pulls the gene annotation
# and terms it is associated to.
#
# param: filename - the gene annotation file
# return: gene_syn - a dictionary; key: a gene, value: array of synonyms
# return: bp_gene_terms - a biological process dictionary; key: a gene,
# value: array of terms that the gene is annotated to
# return: mf_gene_terms - a molecular function dictionary; key: a gene,
# value: array of terms that the gene is annotated to
# return: cc_gene_terms - a cellular component dictionary; key: a gene,
# value: array of terms that the gene is annotated to
def parsing_ann(filename):
    file = open(filename, "r")
    gene_syn = {}
    bp_gene_terms = {}
    mf_gene_terms = {}
    cc_gene_terms = {}

    for line in file:
        if not line.startswith('!'):
            cols = line.strip().split('\t')

            if 'NOT' not in cols[3]:
                gene = cols[2]
                term = cols[4]
                namespace = cols[8]
                synonym_col = cols[10]

                if 'P' in namespace:
                    if gene not in bp_gene_terms.keys():
                        bp_gene_terms[gene] = set()
                    bp_gene_terms[gene].add(term)
                elif 'F' in namespace:
                    if gene not in mf_gene_terms.keys():
                        mf_gene_terms[gene] = set()
                    mf_gene_terms[gene].add(term)
                else:
                    if gene not in cc_gene_terms.keys():
                        cc_gene_terms[gene] = set()
                    cc_gene_terms[gene].add(term)

                if gene not in gene_syn.keys():
                    gene_syn[gene] = set(gene)
                synonyms = synonym_col.split('|')
                for syn in synonyms:
                    if syn not in gene_syn[gene]:
                        gene_syn[gene].add(syn)

    return gene_syn, bp_gene_terms, mf_gene_terms, cc_gene_terms


# This function is for testing the annotation with a modified input file.
#
# param: filename - the file to read in from
# return: genes_terms - dictionary; key: gene, value: terms that the gene is annotated to
def testing_annotation_parsing(filename):
    file = open(filename, "r")
    genes_terms = {}

    for line in file:
        cols = line.split(" ")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]
        if cols[0]:
            genes_terms[cols[0]] = set()
        for i in range(1, len(cols)):
            genes_terms[cols[0]].add(cols[i])

    return genes_terms

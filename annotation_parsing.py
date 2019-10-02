"""
File: annotation_parsing.py
Author: Lily Wise

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
            columns = line.split('\t')
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
            cols = line.split('\t')

            if 'NOT' in cols[3]:
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

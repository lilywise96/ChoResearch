
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
                        bp_gene_terms[gene] = []
                    bp_gene_terms[gene].append(term)
                elif 'F' in namespace:
                    if gene not in mf_gene_terms.keys():
                        mf_gene_terms[gene] = []
                    mf_gene_terms[gene].append(term)
                else:
                    if gene not in cc_gene_terms.keys():
                        cc_gene_terms[gene] = []
                    cc_gene_terms[gene].append(term)

                if gene not in gene_syn.keys():
                    gene_syn[gene] = [gene]
                synonyms = synonym_col.split('|')
                for syn in synonyms:
                    if syn not in gene_syn[gene]:
                        gene_syn[gene].append(syn)

    return gene_syn, bp_gene_terms, mf_gene_terms, cc_gene_terms

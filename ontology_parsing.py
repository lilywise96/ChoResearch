"""
File: ontology_parsing.py
Author: Lily Wise, Joseph Chang

This file parses ontologies for hpo files and for gene ontology files.
"""


# This function parses the hpo ontology file. It pulls the terms and
# their parents to generate a tree.
#
# param: filename - the file that holds the hpo ontology
# return: terms_parents - a dictionary; key: id, value: array of terms (parents)
def hpo_parsing_onto(filename):
    file = open(filename, "r")
    terms_parents = {}
    cur_key = ''

    # Start reading in file.
    for line in file:
        # Find a new term.
        if '[Term]' in line:
            cur_key = ''

        # Find parents.
        elif line.startswith('is_a:'):
            cur_parent = line[6:16]
            if cur_parent not in terms_parents[cur_key]:
                terms_parents[cur_parent] = set()
            terms_parents[cur_key].add(cur_parent)

        # Reads in the id number.
        elif line.startswith('id:'):
            cur_key = line[4:14]
            if cur_key not in terms_parents.keys():
                terms_parents[cur_key] = set()

    return terms_parents


# This function parse the gene ontology file. It pulls the terms and their
# parents. If is_obsolete is found then the term is not included.
#
# param: filename - the file that holds the gene ontology
# return: terms_parents - a dictionary; key: id, value: array of terms (parents)
def parsing_go(filename):
    # TODO: Add part-of terms (currently, only "is-a" terms are added)
    file = open(filename, "r")
    bp_terms_parents = {}
    mf_terms_parents = {}
    cc_terms_parents = {}
    cur_parents = set()
    cur_key = ''
    is_obsolete = False
    namespace = ''

    # Start reading in file.
    end_reached = False
    for line in file:
        line = line.strip()
        # Identifies that a new term is starting.
        if '[Term]' in line:
            if not is_obsolete:
                if namespace is 'b':
                    bp_terms_parents[cur_key] = cur_parents
                elif namespace is 'm':
                    mf_terms_parents[cur_key] = cur_parents
                elif namespace is 'c':
                    cc_terms_parents[cur_key] = cur_parents
            cur_parents = set()
            cur_key = ''
            is_obsolete = False
        # Reads in the id.
        elif line.startswith('id:'):
            cur_key = line[3:].strip()
        # Removes the id if the is_obsolete is found.
        elif 'is_obsolete' in line:
            is_obsolete = True
            if cur_key is not '':
                if namespace is 'b' and cur_key in bp_terms_parents:
                    bp_terms_parents.pop(cur_key)
                elif namespace is 'm' and cur_key in mf_terms_parents:
                    mf_terms_parents.pop(cur_key)
                elif namespace is 'c' and cur_key in cc_terms_parents:
                    cc_terms_parents.pop(cur_key)
        # If it isn't obsolete then the parents can be added if found.
        elif line.startswith('is_a:') and not is_obsolete:
            cur_parents.add(line[6:16])
        elif line.startswith('relationship: part_of') and not is_obsolete:
            cur_parents.add(line[22:32])
        # Checks which namespace it is in.
        elif line.startswith('namespace:'):
            namespace = line[11]
            if namespace is 'b' and cur_key not in bp_terms_parents:
                bp_terms_parents[cur_key] = set()
            elif namespace is 'm' and cur_key not in mf_terms_parents:
                mf_terms_parents[cur_key] = set()
            elif namespace is 'c' and cur_key not in cc_terms_parents:
                cc_terms_parents[cur_key] = set()
        elif line.startswith('[Typedef]') and not end_reached:
            end_reached = True
            if not is_obsolete:
                if namespace is 'b':
                    bp_terms_parents[cur_key] = cur_parents
                elif namespace is 'm':
                    mf_terms_parents[cur_key] = cur_parents
                elif namespace is 'c':
                    cc_terms_parents[cur_key] = cur_parents
            cur_parents = set()
            cur_key = ''
            is_obsolete = False

    # TODO: Resolve and erase debug prints
    print("# of BP terms: " + str(len(bp_terms_parents)))
    print("# of MF terms: " + str(len(mf_terms_parents)))
    print("# of CC terms: " + str(len(cc_terms_parents)))

    # For all terms in an ontology, erase the parents from a different ontology.
    # This step may or may not lose useful information, but it makes further ontology
    # operations less error-prone.
    for term in bp_terms_parents:
        pure_parents = set()
        for parent in bp_terms_parents[term]:
            if parent in bp_terms_parents:
                pure_parents.add(parent)
        bp_terms_parents[term] = pure_parents
    for term in mf_terms_parents:
        pure_parents = set()
        for parent in mf_terms_parents[term]:
            if parent in mf_terms_parents:
                pure_parents.add(parent)
        mf_terms_parents[term] = pure_parents
    for term in cc_terms_parents:
        pure_parents = set()
        for parent in cc_terms_parents[term]:
            if parent in cc_terms_parents:
                pure_parents.add(parent)
        cc_terms_parents[term] = pure_parents

    # GO:0005977 is in BP
    # GO:0006457 is in BP
    # GO:0065009 is in BP
    # They all have KeyError on indexing with tree_child_parent[node] from MF
    
    # Check if parents if each term in an ontology are in the same ontology
    bp_impure_counter = 0
    for term in bp_terms_parents.keys():
        for parent in bp_terms_parents[term]:
            if parent not in bp_terms_parents.keys():
                bp_impure_counter += 1
    mf_impure_counter = 0
    for term in mf_terms_parents.keys():
        for parent in mf_terms_parents[term]:
            if parent not in mf_terms_parents.keys():
                mf_impure_counter += 1
    cc_impure_counter = 0
    for term in cc_terms_parents.keys():
        for parent in cc_terms_parents[term]:
            if parent not in cc_terms_parents.keys():
                cc_impure_counter += 1
    print("# of impure parents in BP: " + str(bp_impure_counter))
    print("# of impure parents in MF: " + str(mf_impure_counter))
    print("# of impure parents in CC: " + str(cc_impure_counter))

    return bp_terms_parents, mf_terms_parents, cc_terms_parents


# Testing of ontology parsing with modified file.
#
# param: filename - the file to parse the ontology from
# return: terms_parents - dictionary; key: term, value: set of parents
def testing_ontology_parsing(filename):
    file = open(filename, "r")
    terms_parents = {}

    for line in file:
        cols = line.split(" ")
        cols[len(cols) - 1] = cols[len(cols) - 1][0:-1]
        if cols[0]:
            terms_parents[cols[0]] = set()
        for i in range(1, len(cols)):
            terms_parents[cols[0]].add(cols[i])

    return terms_parents

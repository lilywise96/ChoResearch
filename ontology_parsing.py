"""
File: ontology_parsing.py
Author: Lily Wise

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
    file = open(filename, "r")
    terms_parents = {}
    cur_parents = set()
    cur_key = ''
    is_obsolete = False

    # Start reading in file.
    for line in file:
        # Identifies that a new term is starting.
        if 'Term' in line:
            if not is_obsolete:
                terms_parents[cur_key] = cur_parents
            cur_parents = set()
            cur_key = ''
            is_obsolete = False
        # Reads in the id.
        elif line.startswith('id:'):
            cur_key = line[4:14]
            if cur_key not in terms_parents.keys():
                terms_parents[cur_key] = set()
        # Removes the id if the is_obsolete is found.
        elif 'is_obsolete' in line:
            is_obsolete = True
            if cur_key is not '':
                terms_parents.pop(cur_key)
        # If it isn't obsolete then the parents can be added if found.
        elif line.startswith('is_a:') and not is_obsolete:
            cur_parents.add(line[6:16])

    return terms_parents


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

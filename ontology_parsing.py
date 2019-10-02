
def hpo_parsing_onto(filename):
    file = open(filename, "r")
    terms_parents = {}
    cur_key = ''

    # Start reading in file.
    for line in file:
        if '[Term]' in line:
            cur_key = ''
        elif line.startswith('is_a:'):
            cur_parent = line[6:16]
            if cur_parent not in terms_parents[cur_key]:
                terms_parents[cur_parent] = []
            terms_parents[cur_key].append(cur_parent)
        elif line.startswith('id:'):
            cur_key = line[4:14]
            if cur_key not in terms_parents.keys():
                terms_parents[cur_key] = []

    return terms_parents


def parsing_go(filename):
    file = open(filename, "r")
    terms_parents = {}
    cur_parents = []
    cur_key = ''
    is_obsolete = False

    for line in file:
        if 'Term' in line:
            if not is_obsolete:
                terms_parents[cur_key] = cur_parents
            cur_parents = []
            cur_key = ''
            is_obsolete = False
        elif line.startswith('id:'):
            cur_key = line[4:14]
            if cur_key not in terms_parents.keys():
                terms_parents[cur_key] = []
        elif 'is_obsolete' in line:
            is_obsolete = True
            if cur_key is not '':
                terms_parents.pop(cur_key)
        elif line.startswith('is_a:') and not is_obsolete:
            cur_parents.append(line[6:16])

    return terms_parents

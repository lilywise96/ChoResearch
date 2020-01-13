"""
Filename: association_creation.py
Author: Lily Wise, Joseph Chang

This creates association and has other functions that are used to calculate the associations.
"""


# This function calculates the coverage given an association and all the itemsets.
#
# param: left_val - the left value of an association
# param: all_itemsets - all of the itemsets given
# return: the number of times the left_val appears in all the itemsets divided by the total number of itemsets
def coverage(left_val, all_itemsets, all_spec):
    count = 0

    # Loop through all the itemsets and check if the left_val is in the itemset
    for itemset in all_itemsets:
        if left_val in all_itemsets[itemset]:
            count += 1

    cover = 0
    if left_val in all_spec:
        cover = count * all_spec[left_val] * 10

    return cover


# This functions calculates the confidence of a given association.
#
# param: all_gt - all of the itemsets
# param: association - the current association to calculate confidence for
#
# returns: the confidence count as a decimal
def confidence(all_gt, association, all_spec):
    # TODO: Find and fix what's causing confidence to be zero for all association rules.
    confidence_count = 0

    # Calculate confidence of an association by iterating over the frequent itemsets

    # Loop through the transactions
    for trans in all_gt:
        found_big = True  # For checking all items in the itemset are present
        # Loop through the items in the itemset
        for associate in association:
            found = False  # For checking just the current item in the itemset is present
            # Loop through each yeast in the transaction
            for i in all_gt[trans]:
                if i == associate:
                    found = True
            if not found:
                found_big = False
        if found_big:
            confidence_count += 1

    conf = 0
    if 1 in association and association[1] in all_spec:
        conf = confidence_count * all_spec[association[1]] * 100

    return conf


# This function takes the frequent itemsets that were created by the apriori algorithm and creates associations. An
# association is kept if it meets the minimum confidence requirements and the left side of the association meets the
# minimum coverage requirements.
#
# param: all_gt - all the itemsets originally read in
# param: freq_itemsets - the frequent itemsets created by the apriori algorithm
# param: min_confidence - the minimum confidence, as a decimal
# param: min_coverage - the minimum coverage, as a decimal
#
# returns: the list of final associations that meets the requirements
def create_associations(left_terms, right_terms, all_gt, freq_itemsets, min_confidence, min_coverage, all_spec):
    print("Starting association creation for min. coverage = " + str(min_coverage)
          + " and min. confidence = " + str(min_confidence))
    print("right terms:")
    print(right_terms)
    final_associations = []
    associations = all_associations(freq_itemsets)

    for associate in associations:
        cur_confidence = confidence(all_gt, associate, all_spec)
        cur_coverage = coverage(associate[0], all_gt, all_spec)

        if cur_confidence >= min_confidence and cur_coverage >= min_coverage \
                and associate[0] in left_terms and associate[1] in right_terms:
            final_associations.append(associate)

    return final_associations


# Creates all possible associations with the frequent itemsets.
#
# param: freq_itemsets - the list of frequent itemsets created by the apriori algorithm
# returns: the associations created.
def all_associations(freq_itemsets):
    associations = []
    for itemset in freq_itemsets:
        association = []
        for item in itemset:
            association.append(item)
        associations.append(association)
        associations.append(association[::-1])

    return associations

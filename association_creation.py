"""
Filename: association_creation.py
Author: Lily Wise


"""
from math import ceil


# This function calculates the coverage given an association and all the itemsets.
#
# param: left_val - the left value of an association
# param: all_itemsets - all of the itemsets given
# return: the number of times the left_val appears in all the itemsets divided by the total number of itemsets
def coverage(left_val, all_itemsets):
    count = 0

    # Loop through all the itemsets and check if the left_val is in the itemset
    for itemset in all_itemsets:
        if left_val in itemset:
            count += 1

    return float(count/len(all_itemsets))


# This functions calculates the confidence of a given association.
#
# param: all_gt - all of the itemsets
# param: association - the current association to calculate confidence for
#
# returns: the confidence count as a decimal
def confidence(all_gt, association):
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

    return float(confidence_count / len(all_gt))


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
def create_associations(all_gt, freq_itemsets, min_confidence, min_coverage):
    final_associations = []
    associations = all_associations(freq_itemsets)
    min_confidence = ceil(min_confidence * len(freq_itemsets))

    for associate in associations:
        cur_confidence = confidence(all_gt, associate)
        cur_coverage = coverage(associate[0], all_gt)
        if cur_confidence > min_confidence and cur_coverage > min_coverage:
            final_associations.append(associate)

    return final_associations


# Creates all possible associations with the frequent itemsets.
#
# param: freq_itemsets - the list of frequent itemsets created by the apriori algorithm
#
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

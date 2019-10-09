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

    for itemset in all_itemsets:
        if left_val in itemset:
            count += 1

    return float(count/len(all_itemsets))


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

    return confidence_count


def all_associations(freq_itemsets):
    associations = []
    for itemset in freq_itemsets:
        association = []
        for item in itemset:
            association.append(item)
        associations.append(association)
        associations.append(association[::-1])

    return associations

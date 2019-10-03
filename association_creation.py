"""
Filename: association_creation.py
Author: Lily Wise


"""
from math import ceil


def create_associations(all_gt, freq_itemsets, min_confidence):
    final_associations = []
    associations = all_associations(freq_itemsets)
    min_confidence = ceil(min_confidence * len(freq_itemsets))

    for associate in associations:
        cur_confidence = confidence(all_gt, associate)
        if cur_confidence > min_confidence:
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

"""
Filename: apriori_algorithm.py
Authors: Joseph Chang, Lily Wise

Calculates the frequent itemsets.
"""

from math import ceil, log10
from itertools import combinations, permutations


# This function calculates support of the itemset from transactions
# param transactions: All transactions in a dictionary
# param itemset: The itemset to calculate support
# return: The support count of the itemset
def weighted_support(transactions, itemset, all_spec):
    support_count = 0

    # Calculate support of an itemset by iterating over the frequent itemsets

    # Loop through the transactions
    for trans in transactions:
        found_big = True  # For checking all items in the itemset are present

        # Loop through the items in the itemset
        for i in itemset:
            found = False  # For checking just the current item in the itemset is present
            # Loop through each yeast in the transaction
            for t in transactions[trans]:
                if i == t:
                    found = True
            if not found:
                found_big = False
        if found_big:
            support_count += 1

    if isinstance(itemset, frozenset):
        itemlist = list(itemset)
        support_weight = 2 * all_spec[itemlist[0]] * all_spec[itemlist[1]] * support_count
        if (all_spec[itemlist[0]] + all_spec[itemlist[1]]) != 0:
            support_weight /= (all_spec[itemlist[0]] + all_spec[itemlist[1]])
        else:
            support_weight = 0
    else:
        support_weight = 2 * all_spec[itemset[0]] * all_spec[itemset[1]] * support_count
        if (all_spec[itemset[0]] + all_spec[itemset[1]]) != 0:
            support_weight /= (all_spec[itemset[0]] + all_spec[itemset[1]])
        else:
            support_weight = 0

    return support_weight


# This function calculates support of the itemset from transactions
# param transactions: All transactions in a dictionary
# param itemset: The itemset to calculate support
# return: The support count of the itemset
def support(transactions, itemset):
    support_count = 0

    # Calculate support of an itemset by iterating over the frequent itemsets

    # Loop through the transactions
    for trans in transactions:
        found_big = True  # For checking all items in the itemset are present

        # Loop through the items in the itemset
        for i in itemset:
            found = False  # For checking just the current item in the itemset is present
            # Loop through each yeast in the transaction
            for t in transactions[trans]:
                if i == t:
                    found = True
            if not found:
                found_big = False
        if found_big:
            support_count += 1

    return support_count


# This function generates a combination from the frequent itemsets of size (itemset_size - 1) and accepts joined
# itemsets if they share (itemset_size - 2) items
# param frequent_itemsets: The table of frequent itemsets discovered
# param itemset_size: The size of joined itemsets
# return: All valid joined itemsets
def generate_selectively_joined_itemsets(frequent_itemsets, itemset_size):

    print("Generating selectively joined itemsets of size " + str(itemset_size) + "...")
    # Record seen_itemsets to prevent duplicates
    seen_itemsets = set()
    joined_itemsets = set()

    # Another implementation...
    # a list of itemsets, eg: [{A, B}, {A, C}, {A, E}, {B, C}, {B, E}, {C, E}]
    fi_list = frequent_itemsets[itemset_size - 1]
    print("size of fi_list: " + str(len(fi_list)))
    # For all pairs of itemsets...
    has_both_ontologies_count = 0
    for i in range(len(fi_list) - 1):
        for j in range(i, len(fi_list)):
            #print(fi_list[i])
            #print(fi_list[j])
            shared_terms = set(fi_list[i]).intersection(set(fi_list[j]))
            #print("shared terms: " + str(shared_terms))
            #print(len(shared_terms))
            # If they share a subset of (size - 2), add the union of the pair to the list.
            if len(shared_terms) >= itemset_size - 2 and len(shared_terms) is not itemset_size - 1:
                #print("itemset 1:")
                #print(fi_list[i])
                #print("itemset 2:")
                #print(fi_list[j])
                joined_set = set(fi_list[i]).union(set(fi_list[j]))
                has_hpo_term = False
                has_go_term = False
                for term in joined_set:
                    if "HP:" in term:
                        has_hpo_term = True
                    if "GO:" in term:
                        has_go_term = True
                if joined_set not in seen_itemsets and has_hpo_term:
                    if has_go_term:
                        has_both_ontologies_count += 1
                    joined_itemsets.add(frozenset(joined_set))
                    seen_itemsets.add(frozenset(joined_set))

    print("For size " + str(itemset_size) + ":")
    print("From selective joining: # of itemsets with both GO and HPO terms: " + str(has_both_ontologies_count))

    # Try all combinations of two itemsets from the table of frequent itemsets and join the pair if they share
    # (itemset_size - 2) items
    # Add each joined itemset to the list if it is not present in the list and discard it otherwise
    # for item1 in frequent_itemsets[itemset_size-1]:
    #     for item2 in frequent_itemsets[itemset_size-1]:
    #
    #         # if the item set is size 1, then you don't need to look for the intersection
    #         if itemset_size-1 == 1:
    #             temp_tuple = (item1, item2)
    #             temp_tuple = tuple(sorted(temp_tuple))
    #             if item1 is not item2 and temp_tuple not in seen_itemsets:
    #                 joined_itemsets.add(temp_tuple)
    #                 seen_itemsets.add(temp_tuple)
    #
    #         # if the item set is greater than 1, then you need to find the intersection
    #         else:
    #             list_a = set(item1)
    #             list_b = set(item2)
    #
    #             # Get the intersection and the union
    #             intersection = list_a.intersection(list_b)
    #             union = list_a.union(list_b)
    #             length_intersection = len(intersection)
    #
    #             # Check if the sets have enough in common
    #             if length_intersection >= itemset_size-2 and length_intersection is not itemset_size-1:
    #                 union = sorted(union)
    #                 temp_tuple = tuple(union)
    #                 has_hpo_term = False
    #                 for term in temp_tuple:
    #                     if "HP:" in term:
    #                         has_hpo_term = True
    #                 if temp_tuple not in seen_itemsets and has_hpo_term:
    #                     seen_itemsets.add(temp_tuple)
    #                     joined_itemsets.add(temp_tuple)

    joined_itemsets = sorted(joined_itemsets)
    return joined_itemsets


# This function checks all the subsets of selected itemsets whether they all
# are frequent or not and prunes the itemset if anyone of the subsets is not frequent
# param selected_itemsets: The itemsets which are needed to be checked
# param frequent_itemsets: The table of frequent itemsets discovered
# param itemset_size: The size of intended frequent itemsets
# return: The itemsets whose all subsets are frequent
def apply_apriori_pruning(selected_itemsets, frequent_itemsets, itemset_size):
    apriori_pruned_itemsets = set()

    # Add each itemset to the list if all of its subsets are frequent and discard it otherwise
    if itemset_size > 3:
        for item in selected_itemsets:
            sub_satisfy = True
            for sub in list(combinations(item, itemset_size-2)):
                if sub not in frequent_itemsets[itemset_size-2]:
                    sub_satisfy = False
            if sub_satisfy:
                apriori_pruned_itemsets.add(item)

    # Add each to the item set if less than 3 because it was already formed from a pruned list so it can't
    # be pruned further.
    else:
        for item in selected_itemsets:
            apriori_pruned_itemsets.add(item)

    # TODO: Remove this debug/print code.
    has_both_ontologies_count = 0
    for itemset in apriori_pruned_itemsets:
        has_hpo_term = False
        has_go_term = False
        for term in itemset:
            if "HP" in term:
                has_hpo_term = True
            if "GO" in term:
                has_go_term = True
        if has_hpo_term and has_go_term:
            has_both_ontologies_count += 1

    print("From apriori pruning: # of itemsets with both GO and HPO terms: " + str(has_both_ontologies_count))

    apriori_pruned_itemsets = sorted(apriori_pruned_itemsets)
    return apriori_pruned_itemsets


# This function generates candidate itemsets of size (itemset_size) by selective joining and apriori pruning
# param frequent_itemsets: The table of frequent itemsets discovered
# param itemset_size: The size of intended frequent itemsets
# return: candidate itemsets formed by selective joining and apriori pruning
def generate_candidate_itemsets(frequent_itemsets, itemset_size):
    joined_itemsets = generate_selectively_joined_itemsets(frequent_itemsets, itemset_size)
    candidate_itemsets = apply_apriori_pruning(joined_itemsets, frequent_itemsets, itemset_size)
    return candidate_itemsets


# This function generates a table of itemsets with all frequent items from transactions based on a given minimum support
# param transactions: The transactions based upon which support is calculated
# param items: The unique set of items present in the transaction
# param min_support: The minimum support to find frequent itemsets
# return: The table of all frequent itemsets of different sizes
def generate_all_frequent_itemsets(transactions, items, min_support, min_weighted_support,
                                   min_information_content, all_spec, all_ic):

    min_support = ceil(min_support * len(transactions))
    min_weighted_support = min_weighted_support / ceil(min_weighted_support * len(transactions))
    min_information_content = min_information_content * -log10(1/len(items))

    frequent_itemsets = dict()
    itemset_size = 0
    frequent_itemsets[itemset_size] = list()
    frequent_itemsets[itemset_size].append(frozenset())

    # Frequent itemsets of size 1
    itemset_size += 1
    frequent_itemsets[itemset_size] = list()

    # Find all frequent itemsets of size 1 and add them to the list
    print("how many items in total: " + str(len(items)))
    count = 0
    test_file = open("test_GO_support.txt", 'w')
    test_file.write("MINIMUM SUPPORT: " + str(min_support) + "\n")
    for item in items:
        #print(str(count))
        count += 1
        list_ver = [item]
        support_check = support(transactions, list_ver)
        # TODO: Check accuracy of support checking
        if "GO" in item:
            test_file.write(item + "\t" + str(support_check))
            test_file.write("\n")
        if support_check >= min_support and all_ic[item] >= min_information_content:
            set_wrapper = set()
            set_wrapper.add(item)
            frequent_itemsets[itemset_size].append(set_wrapper)

    frequent_itemsets[itemset_size] = sorted(frequent_itemsets[itemset_size])

    print("Finished itemsize "+str(itemset_size))
    print(frequent_itemsets[itemset_size])

    # frequent itemsets of greater size
    itemset_size += 1

    # For this study, limit itemset size to 2. One of the items must be an HPO term.
    while frequent_itemsets[itemset_size - 1] and itemset_size <= 2:
        frequent_itemsets[itemset_size] = list()
        candidate_itemsets = generate_candidate_itemsets(frequent_itemsets, itemset_size)
        print("generated candidate itemsets. size = " + str(itemset_size))
        pruned_itemsets = set()

        # Prune the candidate itemset if its support is less than minimum support
        print("size of cand. itemsets: " + str(len(candidate_itemsets)))
        test_count = 0
        for candidate in candidate_itemsets:
            test_count += 1
            #print(test_count)
            #print("candidate info:")
            #print(candidate)
            weighted_sup = weighted_support(transactions, candidate, all_spec)
            #print(weighted_sup)
            if weighted_sup >= min_weighted_support:
                pruned_itemsets.add(candidate)

        print("pruned candidate itemsets for weighted support. size = " + str(itemset_size))
        # TODO: debug/print code.
        has_both_ontologies_count = 0
        for itemset in pruned_itemsets:
            has_hpo_term = False
            has_go_term = False
            for term in itemset:
                if "HP" in term:
                    has_hpo_term = True
                if "GO" in term:
                    has_go_term = True
            if has_hpo_term and has_go_term:
                has_both_ontologies_count += 1
        print("From final support pruning: # of itemsets with both GO and HPO terms: " + str(has_both_ontologies_count))

        pruned_itemsets = list(pruned_itemsets)
        frequent_itemsets[itemset_size] = pruned_itemsets
        print("Finished itemsize " + str(itemset_size))
        print("# of (weighted) frequent itemsets: " + str(len(pruned_itemsets)))
        itemset_size += 1

    return frequent_itemsets


# Calls other methods. The main apriori algorithm.
#
# param: gene_terms - dictionary; key: gene, value: set of terms
# param: gene_set - the set of all distinct genes
# param: min_support - the minimum support
# return: frequent_itemset_table[2] - the frequent itemsets of size 2
def apriori(gene_terms, gene_set, min_support, min_weighted_support, min_information_content,
            all_spec, all_ic):

    frequent_itemset_table = generate_all_frequent_itemsets(gene_terms, gene_set, min_support, min_weighted_support,
                                                            min_information_content, all_spec, all_ic)
    return frequent_itemset_table

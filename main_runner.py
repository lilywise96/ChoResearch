from main import general_main

# def general_main(freq_file_ext, association_file_ext, recreate_onto_ann, recreate_freq_itemsets, tree,
#                  min_support, min_weighted_support, min_confidence, min_information_content, min_coverage):

trees = ['all', 'bp', 'mf', 'hp']
support = [0.02, 0.015]
weighted_support = 0.01
confidence = [0.03, 0.02, 0.01]
info_content = .3
coverage = 0.1

count_freq_file = 1
count_assoc_file = 1

first = True
for tree in range(0, len(trees)):
    for sup in range(0, len(support)):
        for conf in range(0, len(confidence)):
            if first:
                general_main(count_freq_file, count_assoc_file, "true", "true", trees[tree], support[sup],
                             weighted_support, confidence[conf], info_content, coverage)
                first = False
            else:
                general_main(count_freq_file, count_assoc_file, "false", "true", trees[tree], support[sup],
                             weighted_support, confidence[conf], info_content, coverage)
            count_assoc_file += 1
        count_freq_file += 1

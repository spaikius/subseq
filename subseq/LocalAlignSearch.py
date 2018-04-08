"""Description
This module is designed to perform local alignment for each model chain sequence and given target
"""

import SubMatrix
import SmithWaterman


def subseq_la(target, data, matrix, gap_cost, minscore):
    sub_matrix = SubMatrix.SubMatrix(matrix)
    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            aligment = SmithWaterman.SmithWaterman(target
                                                   , sequence
                                                   , gap_cost
                                                   , sub_matrix
                                                   , minscore)

            aligment_list = aligment.get_alignment_coordinates()

            if len(aligment_list) != 0:
                for x, y in aligment_list:
                    start_id = data[model][chain]['ids'][x]
                    end_id = data[model][chain]['ids'][y]

                    match_list.append((model, chain, start_id, end_id))

                print("\n--- Alignment ---")
                print("Model: {0}, chain: {1}".format(model, chain))
                print("Target length: {0} {1}".format(len(target), target[0:40]))
                print("Subject length: {0} {1}".format(len(sequence), sequence[0:40]))

                aligment.print_data()

                print("\n--- END {0} {1} ---".format(model, chain))

    return match_list if len(match_list) != 0 else None

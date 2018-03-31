# @module ss_LA_SEARCH.py
# @public functions: -
# @public variables: -
# @private functions: -
# @private variables: -

import SubMatrix
import SmithWaterman

def subseq_la(target, data, matrix, gap_cost):
    sub_matrix = SubMatrix.SubMatrix(matrix)

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            aligment = SmithWaterman.SmithWaterman(target, sequence, gap_cost, sub_matrix)
            aligment.print_aligment()

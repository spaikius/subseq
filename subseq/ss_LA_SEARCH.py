# @module ss_LA_SEARCH.py
# @public functions: -
# @public variables: -
# @private functions: -
# @private variables: -

import numpy as np


def subseq_la(_target, _data, _matrix, _gap_cost, _extend_cost, _match_value):
    target = _target
    data = _data
    matrix = _matrix
    gap_cost = int(_gap_cost)
    extend_cost = int(_extend_cost)
    match_value = int(_match_value)

    target_header = _create_header_list(target)

    for model in data.keys():
        for chain in data[model].keys():
        	chain_header = _create_header_list(chain)
        	score_matrix = _create_score_matrix(target, chain)

            # _SWA()


def _create_score_matrix(_target, _chain):
	return np.zeros(shape=[len(_target), len(_chain)])


def _create_header_list(_seq):
	return [i for i in _seq.upper()]


def _SWA():
	pass

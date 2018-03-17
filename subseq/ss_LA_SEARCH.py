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

    for model in data.keys():
        for chain in data[model].keys():
            _SWA(target, model, chain, matrix, gap_cost, extend_cost, match_value)


def _SWA():
	pass
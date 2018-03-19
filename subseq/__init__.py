from pymol import cmd

import os
import sys
import re

def __init__(self):
    path = os.path.dirname(__file__)
    sys.path.append(path)

    import ss_ARG
    import ss_DATA
    import ss_RE_SEARCH
    import ss_LA_SEARCH
    import ss_pymol_SELECT


def subseq(*_argv, **_kwargs):
    try:
        ss_ARG.set_parameters(_kwargs)
        param = ss_ARG.get_parameters()

        search_type = param['type']
        target = param['target']
        models = param['models']
        chains = param['chains']
        gap_cost = param['gap']
        extend_cost = param['extend']
        matrix = param['matrix']
        
        data = ss_DATA.get_data(models, chains)

        search_result = None
        if search_type == 're':
            search_result = ss_RE_SEARCH.subseq_re(target, data)
        elif search_type == 'la':
            search_result = ss_LA_SEARCH.subseq_la(target, data, matrix, gap_cost, extend_cost)
        
        if search_result is not None:
            ss_pymol_SELECT.select(search_result)
        else:
            print("Empty")

    except Exception as e:
        print(' Error: ' + str(e))

    ss_ARG.reset_parameters()



cmd.extend('subseq', subseq)

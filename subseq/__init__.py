from pymol import cmd

import os
import sys
import re

def __init__(self):
    path = os.path.dirname(__file__)
    sys.path.append(path)

    import ParseKwargs
    import Data
    import RegExSearch
    import LocalAlignSearch
    import Select


def subseq(*_argv, **_kwargs):
    try:
        ParseKwargs.set_parameters(_kwargs)
        param = ParseKwargs.get_parameters()

        search_type = param['type']
        target = param['target']
        models = param['models']
        chains = param['chains']
        gap_cost = param['gap']
        extend_cost = param['extend']
        matrix = param['matrix']
        
        data = Data.get_data(models, chains)

        search_result = None
        if search_type == 're':
            search_result = RegExSearch.subseq_re(target, data)
        elif search_type == 'la':
            search_result = LocalAlignSearch.subseq_la(target, data, matrix, gap_cost, extend_cost)
        
        if search_result is not None:
            Select.select(search_result)
        else:
            print("Empty")

    except Exception as e:
        print(' Error: ' + str(e))

    ParseKwargs.reset_parameters()



cmd.extend('subseq', subseq)

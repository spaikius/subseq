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
        parameters = ParseKwargs.ParseKwargs(_kwargs)

        data = Data.Data(parameters['models'], parameters['chains'])

        search_result = None
        if parameters['algorithm'] == 're':
            search_result = RegExSearch.subseq_re(parameters['target'], data)
        elif parameters['algorithm'] == 'la':
            search_result = LocalAlignSearch.subseq_la(parameters['algorithm'], data, parameters['submatrix'], parameters['gapcost'], '2')
        
        if search_result is not None:
            Select.select(search_result)
        else:
            print("Empty")
    except Exception as e:
        print(' Error: ' + str(e))


cmd.extend('subseq', subseq)

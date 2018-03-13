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
    import ss_pymol_SELECT


def subseq(*_argv, **_kwargs):
    try:
        ss_ARG.set_parameters(_kwargs)
        param = ss_ARG.get_parameters()

        data = ss_DATA.get_data(param['models'], param['chains'])

        search_result = None
        if param['type'] == 're':
            search_result = ss_RE_SEARCH.subseq_re(param['target'], data)

        ss_pymol_SELECT.select(search_result)
    except Exception as e:
        print(' Error: ' + str(e))

    ss_ARG.reset_parameters()



cmd.extend('subseq', subseq)


subseq_help = '''
@usage : subseq target=DIEVDLLKNGER, type=gl, chains=[A,C,D]
@@ Please note: arguments must be sepereted with comma (,) 
@params:
(required) target=(str)   : target sequence
(optional) model=[array]  : default all
(optional) type=(str)     : re or lo default re
(optional) chains=[array] : default all
'''

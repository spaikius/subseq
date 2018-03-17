# @module ss_pymol_select.py
# @public functions: select
# @public variables: -
# @private functions: -
# @private variables: -

from pymol import cmd
from random import randint


# @function type: public
# @arguments: list from ss_RE_SEARCH.subseq_re() or ss_LA_SEARCH.subseq_la()
# @returns: -
# @description: for given data creates a 'ss_%{randint}' selection and fills it with residues
def select(_select_list):
    _select_list = _select_list
    lower_bound = 0
    upper_bound = 10000
    # select ID
    select_name = 'ss_' + str(randint(lower_bound, upper_bound))
    # empty select
    cmd.select(select_name, None) 

    for sele_tuple in _select_list:
        model = sele_tuple[0]
        chain = sele_tuple[1]
        start_id = sele_tuple[2]
        end_id = sele_tuple[3]
        
        # select /model/?/chain/resi
        cmd.select(select_name, select_name +  
            " | /" + str(model) +  
            "//" + str(chain) + 
            "/" + str(start_id) + 
            "-" + str(end_id))

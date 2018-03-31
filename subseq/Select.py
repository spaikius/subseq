"""Description

This module is designed to execute pymol select querry

Module attribute:
    select(select_list: list)
        workflow:
            1) create a random select id
            2) create an empty select
            3) for each tuple in select_id:
                a) construct select querry from data in select_list
                   (function parameter)
                b) execute and append selection to select_id

"""

from pymol import cmd
from random import randint


def select(select_list):
    lower_bound = 0
    upper_bound = 100000
    # random select ID
    select_id = 'ss_' + str(randint(lower_bound, upper_bound))
    # empty select
    cmd.select(select_id, None)

    select_query = None

    for sele_tuple in select_list:
        model = sele_tuple[0]
        chain = sele_tuple[1]
        start_id = sele_tuple[2]
        end_id = sele_tuple[3]

        # select /model/?/chain/resi
        select_query = " | /{0}//{1}/{2}-{3}".format(model, chain, start_id, end_id)

        # Execute and append selection to select_id
        cmd.select(select_id, select_id + select_query)

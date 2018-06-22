import string
import re

from pymol import cmd, stored

stored.id = 0

def select(select_list, target, sele, method):
    """Creates pymol selection object"""

    select_name = string.Formatter().vformat(
        sele,
        (),
        SafeDict(
            method=method,
            target=re.sub(r'[^\w]', '', target)[:10],
            id=new_id(sele)))

    # empty select
    cmd.select(select_name, None)

    for select_tuple in select_list:
        model = select_tuple[0]
        chain = select_tuple[1]
        resi = select_tuple[2]

        # select /model/?/chain/resi-resi
        select_query = " | /{0}//{1}/{2}".format(model, chain, resi)

        # Execute and append selection to select_id
        cmd.select(select_name, select_name + select_query)


def new_id(sele):
    """Returns an incremented id if id token is requested"""
    if re.search(r'{id}', sele):
        stored.id += 1

    return stored.id


class SafeDict(dict):
    def __missing__(self, key):
        return key

"""Description
This module is designed to search requested subsequence using Regular Expresion
This module provides a function for searching 
"""
import re
from Exceptions import BadRegExSyntaxError


def subseq_re(target, data):
    """
    workflow:
        1) create a RegExp object from target (function parameter)
        2) scan data (function parameter) by using RegExp object
        3) append information about match to match_list
        4) return match_list if its length is not 0 else return None
    """
    sequence = 'sequence'
    ids = 'ids'
    match_list = list()

    # RegExp validation
    try:
        # re.I - ignore case sensetive
        re_target = re.compile(target, re.I)
    except:
        raise BadRegExSyntaxError('Bad syntax - target(RegExp): ' + str(target))

    # scan data by using RegExp object
    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain][sequence]):

                start_id = data[model][chain][ids][match.start()]
                end_id = data[model][chain][ids][match.end() - 1]
                
                match_list.append((model, chain, start_id, end_id))

    # return match list if its length is not 0 else return None
    return match_list if len(match_list) != 0 else None

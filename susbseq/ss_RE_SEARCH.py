# @module ss_RE_SEARCH.py
# @public functions: subseq_re
# @public variables: -
# @private functions: -
# @private variables: -

import re


# @function type: public
# @arguments: target, data (from ss_DATA.get_data())
# @returns: list of tuple (model, chain, start_pos, end_pos)
#           or 'None' if no matches have been found
# @description: for given target seacrh through data using RegExp method
# warn if no matche have been found
def subseq_re(_target, _data):
    target = _target
    data = _data
    sequence = 'sequence'
    ids = 'ids'
    match_list = list()

    # RegExp validation
    try:
        # re.I - ignore case sensetive
        re_target = re.compile(target, re.I)
    except:
        raise Exception('Bad syntax - target(RegExp): ' + str(_target))

    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain][sequence]):
                start_id = data[model][chain][ids][match.start()]
                end_id = data[model][chain][ids][match.end() - 1]
                match_list.append((model, chain, start_id, end_id))

    if len(match_list) != 0:
        return match_list
    else:
        print("Warrning: No mathces have been found for given target: " + str(target))
        return None
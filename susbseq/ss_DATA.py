# @module ss_DATA.py
# @public functions: get_data
# @public variables: -
# @private functions: _init_data,
#                     _replace_to_one_letter,
#                     _create_empty_dict,
#                     _fill_dict,
#                     _get_chains
# @private variables: _one_letter

from pymol import cmd

_one_letter = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'CA': ''
}


# @function type: public
# @arguments: model list
# @returns: nested dictionary {model: {chain: sequence, },  }
# @description: encapsulation, returns nested dictionary
def get_data(_models, _chains):
    # get all data from pymol
    aa_dict = _init_data()
    # replace all 3 letters
    _replace_to_one_letter(aa_dict) 
     # creates an empty dict for complete_chains
    complete_chains = _create_empty_dict(_models, _chains)
    # using pymol data fills complete_chains with data
    _fill_dict(_models, _chains, complete_chains, aa_dict) 

    return complete_chains


# @function type: private
# @arguments: -
# @returns: dict
# @description: return aa_dict{[...]} with information about residue 
def _init_data():
    aa_dict = dict()
    aa_dict['aa_list'] = list()

    # iterate through all c-alpha atoms and return
    # position number, three letter code, corresponding chain and model
    cmd.iterate("(name ca)",
                "aa_list.append([resn, resi, chain, model])",
                space=aa_dict)
    return aa_dict


# @function type: private
# @arguments: -
# @returns: -
# @description: replaces all 3 letters code aa to single letter code aa
def _replace_to_one_letter(_aa_dict):
    for aa_list in _aa_dict['aa_list']:
        aa_list[0] = _one_letter[aa_list[0]]


# @function type: private
# @arguments: -
# @returns: -
# @description: constructs an empty nested dictionary _complete_chains
# {model: {chain: {sequence: str, ids: list, }, },  }
def _create_empty_dict(_models, _chains):
    complete_chains = dict()
    sequence = 'sequence'
    ids = 'ids'

    for model in _models:
        complete_chains[model] = dict()

        for chain in _get_chains(model):
            # skip if chain is not requested
            if chain not in _chains:
                continue

            complete_chains[model][chain] = dict()
            complete_chains[model][chain][sequence] = ''
            complete_chains[model][chain][ids] = list()

    return complete_chains


# @function type: private
# @arguments: -
# @returns: -
# @description: constructs a nested dictionary _complete_chains 
# {model: {chain: {sequence: str, ids: list, }, },  }
# warns if any empty model is presented (does not contain any chain from users requested chains)
def _fill_dict(_models, _chains, _complete_chains, _aa_dict):
    sequence = 'sequence'
    ids = 'ids'

    for resn, resi, chain, model in _aa_dict['aa_list']:
        # if model or chain is not requested, skip iteration
        if model not in _models or chain not in _chains:
            continue

        _complete_chains[model][chain][sequence] += resn
        _complete_chains[model][chain][ids].append(resi)

    # remove all empty models and warn if found any 
    for model in _complete_chains.keys():
        if not _complete_chains[model]:
            print(" Warning: model '" + str(model) + "' does not contain any requested chain\n" +
                "  Models '" + str(model) + "' chains: " + str(_get_chains(model)) + '\n' +
                "  Requested chains: " + str(_chains))
            del _complete_chains[model]


# @function type: private
# @arguments: a list of models
# @returns: a list of chains with no duplicates
# @description: returns a list of chains in model
def _get_chains(_model):
    chain_list = list()
    for chain in cmd.get_chains(_model):
            chain_list.append(chain)

    return chain_list

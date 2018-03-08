# @module ss_DATA.py
# @public functions: get_data,
#                    del_data
# @public variables: -
# @private functions: _init_data,
#                     _filter_data,
#                     _construct_aa_chains
# @private variables: _aa_dict, _complete_chains

from pymol import cmd

_oneLetter = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'CA': ''
}

_aa_dict = dict()
_complete_chains = dict()


# @function type: public
# @arguments: model list, chain list
# @returns: nested dictionary {model: {chain: sequence, },  }
# @description: encapsulation, returns nested dictionary
def get_data(_models, _chains):
    _init_data()
    _filter_data(_models, _chains)
    _construct_aa_chains(_models, _chains)

    return _complete_chains


# @function type: public
# @arguments: -
# @returns: -
# @description: frees memory
def del_data():
    try:
        del _aa_dict
        del _complete_chains
    except NameError:
        return


# @function type: private
# @arguments: -
# @returns: -
# @description: fills _aa_dict{[...]} with information about residue 
def _init_data():
    _aa_dict['aa_list'] = list()
    # iterets through all c-alpha atoms and returns
    # possition number, three letter code, coresponding chain and model
    cmd.iterate("name ca",
                "aa_list.append([resn, resi, chain, model])",
                space=_aa_dict)


# @function type: private
# @arguments: model list, chain list
# @returns: -
# @description: filters _aa_dict 
def _filter_data(_models, _chains):
    aa_list = _aa_dict['aa_list']

    for i in range(0, len(aa_list)):
        if aa_list[i][3] not in _models or \
                        aa_list[i][2] not in _chains:
            # delete list element
            del _aa_dict['aa_list'][i]


# @function type: private
# @arguments: model list, chain list
# @returns: -
# @description: constructs a nested dictionary _complete_chains 
# {model: {chain: sequence, },  }
def _construct_aa_chains(_models, _chains):
    for aa in _aa_dict['aa_list']:
        aa[0] = _oneLetter[aa[0]]

    for model in _models:
        _complete_chains[model] = dict()
        for chain in _chains:
            _complete_chains[model][chain] = ''

    for i in _aa_dict['aa_list']:
        _complete_chains[i[3]][i[2]] += str(i[0])

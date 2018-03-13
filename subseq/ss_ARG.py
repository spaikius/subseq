# @module ss_ARG.py
# @public functions: set_parameters,
#                    get_paramaters,
#                    reset_parameters
# @public variables: -
# @private functions: _set_parameters,
#                     _validate_parameters,
#                     _parse_str_to_list,
#                     _get_all_modules,
#                     _get_all_chains
# @private variables: _algorithm_type_values, _params

from pymol import cmd
import re

_algorithm_type_values = ['re',  # Regular Expresion (default)
                'la'   # Local alignment
                ]

_params = {
    'type': None,
    'target': None,
    'models': None,
    'chains': None,
}


# @function type: public
# @arguments: pointer to a dictionary
# @returns: -
# @description: encapsulation
def set_parameters(_kwargs):
    _set_parameters(_kwargs)
    _validate_parameters()


# @function type: public
# @arguments: -
# @returns: dictionary
def get_parameters():
    return _params


# @function type: public
# @arguments: -
# @returns: -
# @description: sets all values of a dictionary (_params) to 'None'
def reset_parameters():
    for key in _params.keys():
        _params[key] = None


# @functions type: private
# @arguments: pointer to a dictionary
# @returns: -
# @description: sets all values of a dictionary (_params)
# to a coresponding given key
def _set_parameters(_kwargs):
    # Remove '_self' key if presented
    _kwargs.pop('_self', None)
    for key, value in _kwargs.items():
        key = key.lower()

        # Check if given keys exists in _params dictionary
        # if not raise an exception
        if key in _params:
            _params[key] = value
        else:
            raise Exception('Unknown parameter: ' + str(key))


# @functions type: private
# @arguments: -
# @returns: -
# @description: validates all values of a dictionary (_params).
# If key value is not provided, it will set the default value
def _validate_parameters():
    algorithm_type = _params['type']
    default_algo_type = _algorithm_type_values[0]
    target = _params['target']
    models = _params['models']
    chains = _params['chains']

    pymol_models = _get_all_models()
    pymol_models_chains = _get_all_chains(pymol_models)

    #-- Algorithm type validation
    if algorithm_type is None:
        _params['type'] = default_algo_type
    else:
        if algorithm_type not in _algorithm_type_values:
            raise Exception("Undefined algorithm type: " + str(algorithm_type))

    #-- Target validation
    if target is None:
        raise Exception("Parameter 'target=' must be defined")

    #-- Models validation
    if models is None:
        _params['models'] = pymol_models
    else:
        models_list = _parse_str_to_list(models)
        for model in models_list:
            if model not in pymol_models:
                raise Exception("Model '" + str(model) + "' is not available\n" +
                    "Available models: " + str(pymol_models))

        _params['models'] = models_list

    #-- Chains validation
    if chains is None:
        _params['chains'] = pymol_models_chains
    else:
        chains_list = _parse_str_to_list(chains)
        for chain in chains_list:
            if chain not in pymol_models_chains:
                raise Exception("Chain '" + str(chain) + "' is not available\n" +
                    "Available chains: " + str(pymol_models_chains))

        _params['chains'] = chains_list

# @function type: private
# @arguments: string
# @returns: list
# @descriptions: parses string from _kwargs values and returns list 
def _parse_str_to_list(_str):
    # split by any seperator (, . / etc.)
    raw_str = re.sub('([A-Za-z0-9_]+)', r'\1', _str)
    # remove symbols '[', ']', ''', '"', white spaces
    raw_str = re.sub('[\[\]\'\"\s]', '', raw_str)
    return raw_str.split(',')


# @function type: private
# @arguments: -
# @returns: a list of models
# @descriptions: returns names list of all existing models in pymol
# raises an exception if no model was found
def _get_all_models():
    models = cmd.get_names()

    if not bool(models):
        raise Exception("No models are currently opened")

    return models
  

# @function type: private
# @arguments: a list of models
# @returns: a list of chains with no duplicates
# @description: gets all chains in provided models and returns a list
# without duplicates
def _get_all_chains(_models):
    chain_list = list()
    for model in _models:
        for chain in cmd.get_chains(model):
            chain_list.append(chain)

    return list(set(chain_list))

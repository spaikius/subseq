# @module ss_ARG.py
# @public functions: set_parameters,
#                    get_paramaters,
#                    reset_parameters
# @public variables: -
# @private functions: _set_parameters,
#                     _validate_parameters,
#                     _get_modules
#                     _get_chains
# @private variables: _type_values, _params

from pymol import cmd

_type_values = ['re',  # Regular Expresion (default)
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
            raise Exception('Unknown parameter: ', key)


# @functions type: private
# @arguments: -
# @returns: -
# @description: validates all values of a dictionary (_params).
# If key value is not provided, it will set the default value
def _validate_parameters():
    algorithm_type = _params['type']
    target = _params['target']
    models = _params['models']
    chains = _params['chains']

    pymol_models = _get_models()
    pymol_models_chains = _get_chains(pymol_models)

    if algorithm_type is None:
        # default value: re
        _params['type'] = _type_values[0]
    else:
        if algorithm_type not in _type_values:
            raise Exception("Unknown 'type' value: ", _params['type'])

    # --
    if target is None:
        raise Exception("Parameter 'target' must be defined")
    elif not isinstance(target, str):
        raise Exception("Parameter 'target' bust be type of 'str'")

    # --
    if models is None:
        # default value: all existing models
        _params['models'] = _get_models()
    # only 1 model is given
    elif isinstance(models, str):
        # check if it exists
        if models not in pymol_models:
            raise Exception('Model: ', models, ' is undefined')
    # a list of models is given
    elif isinstance(models, list):
        for model in models:
            # check if model type is str
            if not isinstance(model, str):
                raise Exception('Models value can not be a nested list')
            # check if model exists
            if model not in pymol_models:
                raise Exception('Models: ', models, ' are undefined')
    # the given object is not a type of str or list
    else:
        raise Exception("models value must be type of 'str' or 'list'")

    # --
    if chains is None:
        # default value: all existing chains
        _params['chains'] = pymol_models_chains
    else:
        if isinstance(chains, list):
            required_chains = _get_chains(models)
            for chain in chains:
                if chain not in required_chains:
                    raise Exception('Chain: ', chain, 'can not be',
                                    ' found in: ', models)
        else:
            raise Exception("chains value must be type of 'list'")


# @function type: private
# @arguments: -
# @returns: a list of models
# @descriptions: returns names list of all existing models in pymol
# raises an exception if no model was found
def _get_models():
    models = cmd.get_names()

    if not bool(models):
        raise Exception("No models are currently opened")

    return models
  

# @function type: private
# @arguments: a list of models
# @returns: a list of chains with no duplicates
# @description: gets all chains in provided models and returns a list
# without duplicates
def _get_chains(_models):
    chain_list = list()
    for model in _models:
        for chain in cmd.get_chains(model):
            chain_list.append(chain)

    return list(set(chain_list))

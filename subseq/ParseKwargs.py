"""Description

This module is designed to parse and validate raw kwargs (keyword arguments)
from Pymol command line.
    
ParseKwargs (class):
    class description:
        This class is designed to parse and validate
        raw kwargs (keyword arguments)    

    class param attribute keys:
        algorithm, target, models, chains, gapcost, submatrix, minscore

    class attributes:
        - algorithm_types (list): static.
            All avialable algorithms 

        - default_algorithm (str): static.
            Default algorithm (re) RegExp

        - default_gap_cost (int): static.
            Default gap cost (10)

        - default_score_matrix (str): static.
            Default score matrix (blossum62)

        - default_minscore (float): static.
            Default minimum score for local aligment (51.0%)

        - parameters (dict): typed.
            Contains all parameters values
        
        - __init__ (kwargs: dict): typed.
            Constructor

        - __getitem__ (param: str) -> str: typed.
            Getter. Returns a value for given key from `param` dictionary
                
        - set_parameters (kwargs: dict): typed.
            Set parameters(class attribute) key-value pairs

        - validate_all_parameters (): typed.
            Calls all *_validation() functions

        - algorithm_validation (): typed.
            Validates algorithm type

        - target_validation (): typed.
            Validates target

        - models_validation (): typed.
            Validates models

        - chains_validation (): typed.
            Validates chains

        - gapcost_validation (): typed.
            Validates gap cost

        - submatrix_validation (): typed.
            Validates substitution matrix

        - minscore_validation (): typed.
            Validates minimum score (in percentages) for local aligment

Example:
    ParseKwargs_obj = ParseKwargs(kwargs)
    target_value = ParseKwargs_obj['target']

"""

import os

from HelperFunctions import get_all_models, get_all_chains, parse_str_to_list

class ParseKwargs:
    algorithm_types = [
        're',  # Regular Expression (default)
        'la'   # Local alignment
    ]

    default_algorithm = algorithm_types[0]
    default_gap_cost = 10
    default_score_matrix = 'blosum62'
    default_minscore = 51.

    def __init__(self, kwargs):
        self.parameters = {
            'algorithm':   None,
            'target':      None,
            'models':      None,
            'chains':      None,
            'gapcost':     None,
            'submatrix':   None,
            'minscore':    None,
        }

        self.set_parameters(kwargs)
        self.validate_all_parameters()

    def __getitem__(self, param):
        return self.parameters[param]

    def set_parameters(self, kwargs):
        # Remove '_self' key if presented
        kwargs.pop('_self', None)

        # Set self.parameters[key] values from kwargs key-value pairs
        for key, value in kwargs.items():
            key = key.lower()

            if key in self.parameters:
                self.parameters[key] = value
            else:
                # If there is no such key in self.parameters, then raise exception
                raise Exception("Unexpected parameter: {}".format(key))

    def validate_all_parameters(self):
        self.algorithm_validation()
        self.target_validation()
        self.models_validation()
        self.chains_validation()
        self.gapcost_validation()
        self.submatrix_validation()
        self.minscore_validation()

    def algorithm_validation(self):
        if self.parameters['algorithm'] is None:
            self.parameters['algorithm'] = self.default_algorithm
        else:
            if self.parameters['algorithm'] not in self.algorithm_types:
                raise Exception("Unexpected algorithm type {}"
                                .format(self.parameters['algorithm']))

    def target_validation(self):
        if self.parameters['target'] is None:
            raise Exception("Parameter 'target=' must be defined"
                            .format(self.errors['targetValueError']))
        else:
            self.parameters['target'] = self.parameters['target'].upper()

    def models_validation(self):
        all_pymol_models = get_all_models()

        if self.parameters['models'] is None:
            self.parameters['models'] = all_pymol_models
        else:
            models_list = parse_str_to_list(self.parameters['models'])
            models_list = [model.upper() for model in models_list]

            # Check if all models are opened by pymol
            for model in models_list:
                if model not in all_pymol_models:
                    raise Exception("Found unknown model {}".format(model))

            self.parameters['models'] = models_list   

    def chains_validation(self):
        all_models_chains = get_all_chains(self.parameters['models'])

        if self.parameters['chains'] is None:
            self.parameters['chains'] = all_models_chains
        else:
            chains_list = parse_str_to_list(self.parameters['chains'])
            chains_list = [chain.upper() for chain in chains_list]
            for chain in chains_list:
                if chain not in all_models_chains:
                    raise Exception("Found unknown chain {}".format(chain))

            self.parameters['chains'] = chains_list

    def gapcost_validation(self):
        if self.parameters['gapcost'] is None:
            self.parameters['gapcost'] = self.default_gap_cost
        else:
            try:
                gap_cost = float(self.parameters['gapcost'])
                self.parameters['gapcost'] = gap_cost
            except:
                raise Exception("Gap cost must be \
                                type of int or float. Got {}"
                                .format(self.parameters['gapcost']))

    def submatrix_validation(self):
        if self.parameters['submatrix'] is None:
            self.parameters['submatrix'] = \
                os.path.join(os.path.dirname(__file__)
                             , self.default_score_matrix)
        else:
            if not os.path.isfile(self.parameters['submatrix']):
                raise Exception("File or dir does not exist {}"
                                .format(os.path.basename(
                                    self.parameters['submatrix'])))

    def minscore_validation(self):
        if self.parameters['minscore'] is None:
            self.parameters['minscore'] = self.default_minscore
        else:
            try:
                minscore = float(self.parameters['minscore'])
                if minscore >= 0 and minscore <= 100:
                    self.parameters['minscore'] = minscore
                else:
                    raise Exception('minscore value must be between 0 and 100')

            except:
                raise Exception('minscore should be type of int or float')

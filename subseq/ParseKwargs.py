"""Description

This module provides a class for parsing and validating raw kwargs (keyword arguments)
from Pymol command line.
    
"""

import os

from HelperFunctions import get_all_models, get_all_chains, parse_str_to_list
from Exceptions import BadParameterError


class ParseKwargs:
    method_types = [
        're',  # Regular Expression (default)
        'la'   # Local alignment
    ]

    default_method = method_types[0]
    default_gap_cost = 10
    default_score_matrix = os.path.join(os.path.dirname(__file__), 'blosum62')
    default_minscore = 51.

    def __init__(self, kwargs):
        self.parameters = {
            'method':   None,
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
                # raise exception If there is no such key in self.parameter
                raise BadParameterError("Unexpected parameter: {}".format(key))

    def validate_all_parameters(self):
        self.method_validation()
        self.target_validation()
        self.models_validation()
        self.chains_validation()
        self.gapcost_validation()
        self.submatrix_validation()
        self.minscore_validation()

    def method_validation(self):
        # If method is not provided set its value to default (Regular Expression)
        if self.parameters['method'] is None:
            self.parameters['method'] = self.default_method
        else:
            # check if requested algortihmn exists
            if self.parameters['method'] not in self.method_types:
                raise BadParameterError("Unexpected method type: {}"
                                .format(self.parameters['method']))

    def target_validation(self):
        # Raise exception if target is not provided
        if self.parameters['target'] is None:
            raise BadParameterError("Parameter 'target=' must be defined")
        else:
            self.parameters['target'] = self.parameters['target'].upper()

    def models_validation(self):
        all_pymol_models = get_all_models()

        if self.parameters['models'] is None:
            self.parameters['models'] = all_pymol_models
        else:
            models_list = parse_str_to_list(self.parameters['models'])
            models_list = [model.upper() for model in models_list]

            # Check if all models exists
            for model in models_list:
                if model not in all_pymol_models:
                    raise BadParameterError("Found unknown model: {}".format(model))

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
                    raise BadParameterError("Found unknown chain: {}".format(chain))

            self.parameters['chains'] = chains_list

    def gapcost_validation(self):
        if self.parameters['gapcost'] is None:
            self.parameters['gapcost'] = self.default_gap_cost
        else:
            try:
                gap_cost = float(self.parameters['gapcost'])
            except:
                raise BadParameterError("Gap cost must be"
                                + "type of int or float. Got {}"
                                  .format(self.parameters['gapcost']))

            if gap_cost < 0.1:
                raise BadParameterError("Gap cost can not be less than 0.1")
            self.parameters['gapcost'] = gap_cost

    def submatrix_validation(self):
        if self.parameters['submatrix'] is None:
            self.parameters['submatrix'] = self.default_score_matrix
        else:
            # Check if file exists.
            if not os.path.isfile(self.parameters['submatrix']):
                raise BadParameterError("File or dir does not exist: {}"
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
                    raise BadParameterError('minscore value must be between 0 and 100')

            except:
                raise BadParameterError('minscore should be type of int or float')

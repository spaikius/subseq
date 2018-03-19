"""Description

This module is designed to parse and validate raw kwargs (keyword arguments)
from Pymol command line.

"""

import os

from HelperFunctions import get_all_models, get_all_chains, parse_str_to_list

class ParseKwargs:
    errors = {
        'inputError':         '1',
        'algorithmTypeError': '2',
        'targetValueError':   '3',
        'unknownModelError':  '4',
        'unknownChainError':  '5',
        'gapCostError':       '6',
        'fileError':          '7',
    }

    algorithm_types = [
        're',  # Regular Expression (default)
        'la'   # Local alignment
    ]

    default_algorithm = algorithm_types[0]
    default_gap_cost = 10
    default_score_matrix = 'blosum62'

    def __init__(self, kwargs):
        self.parameters = {
            'algorithm':   None,
            'target':      None,
            'models':      None,
            'chains':      None,
            'gapcost':     None,
            'submatrix':   None,
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
            value = value.lower()

            if key in self.parameters:
                self.parameters[key] = value
            else:
                # If there is no such key in self.parameters, then raise exception
                raise Exception("Error number: {}, unexpected parameter: {}"
                                .format(self.errors['inputError'], key))

    def validate_all_parameters(self):
        self.algorithm_validation()
        self.target_validation()
        self.models_validation()
        self.chains_validation()
        self.gapcost_validation()
        self.submatrix_validation()

    def algorithm_validation(self):
        if self.parameters['algorithm'] is None:
            self.parameters['algorithm'] = self.default_algorithm
        else:
            if self.parameters['algorithm'] not in self.algorithm_types:
                raise Exception("Error number: {}, unexpected algorithm type {}"
                                .format(self.errors['algorithmTypeError']
                                        , self.parameters['algorithm']))

    def target_validation(self):
        if self.parameters['target'] is None:
            raise Exception("Error number: {}, parameter 'target=' must be defined"
                            .format(self.errors['targetValueError']))
        else:
            self.parameters['target'] = self.parameters['target'].upper()

    def models_validation(self):
        all_pymol_models = get_all_models()

        if self.parameters['models'] is None:
            self.parameters['models'] = all_pymol_models
        else:
            models_list = parse_str_to_list(self.parameters['models'])

            # Check if all models are opened by pymol
            for model in models_list:
                if model not in all_pymol_models:
                    raise Exception("Error number: {}, found unknown model {}"
                                    .format(self.errors['unknownModelError']
                                            , model))

            self.parameters['models'] = models_list   

    def chains_validation(self):
        all_models_chains = get_all_chains(self.parameters['models'])

        if self.parameters['chains'] is None:
            self.parameters['chains'] = all_models_chains
        else:
            chains_list = parse_str_to_list(self.parameters['chains'])
            for chain in chains_list:
                if chain not in all_models_chains:
                    raise Exception("Error number: {}, found unknown chain {}"
                                    .format(self.errors['unknownChainError']
                                            , chain))

            self.parameters['chains'] = chains_list

    def gapcost_validation(self):
        if self.parameters['gapcost'] is None:
            self.parameters['gapcost'] = self.default_gap_cost
        else:
            try:
                gap_cost = float(self.parameters['gapcost'])
                self.parameters['gapcost'] = gap_cost
            except:
                raise Exception("Error number: {}, gap cost must be \
                                type of int or float. Got {}"
                                .format(self.errors['gapCostError']
                                        , self.parameters['gapcost']))

    def submatrix_validation(self):
        if self.parameters['submatrix'] is None:
            self.parameters['submatrix'] = \
                os.path.join(os.path.dirname(__file__)
                             , self.default_score_matrix)
        else:
            if not os.path.isfile(self.parameters['submatrix']):
                raise Exception("Error number: {}, file or dir does not exist {}"
                                .format(self.errors['fileError']
                                        ,os.path.basename(self.parameters['submatrix'])))

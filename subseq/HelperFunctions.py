"""Description

This module contains only helper functions:
get_all_models() -> list 				; return all pymol models
get_all_chains(models: list) -> list 	; return all chains for given models list
get_chains(model: str) -> list 			; return all chains for given model
parse_str_to_list(string: str) -> list 	; parse a string and return a list of values

"""

import re
from pymol import cmd

def get_all_models():
    """ return all pymol models """
    
    models = [model.upper() for model in cmd.get_names()]

    # if a list is empty -> no models are currently avialable
    if not bool(models):
        raise Exception("No models are currently opened")

    return models
  

def get_all_chains(models):
    """ return a list of chains for given models """

    chain_list = list()

    for model in models:
        for chain in cmd.get_chains(model):
            chain_list.append(chain.upper())

    return list(set(chain_list))


def get_chains(model):
    """ return a list of chains for given only one model """

    chain_list = list()
    for chain in cmd.get_chains(model):
            chain_list.append(chain.upper())

    return chain_list


def parse_str_to_list(string):
    """ Parse a string and return a list of values """

    # split by any seperator (, . / etc.) excluding all letters, numbers and _
    raw_str = re.sub('([A-Za-z0-9_]+)', r'\1', string)

    # remove symbols '[', ']', ''', '"', white space
    raw_str = re.sub('[\[\]\'\"\s]', '', raw_str)

    return raw_str.split(',')

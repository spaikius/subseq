import logging
import re
import os

from pymol import cmd

def parse_targets(targets):
    """Parser for user input"""
    if targets is '':
        logging.error("parameter 'targets' is not specified.")
        return

    targets_list = list()

    for target in targets.split(" "):

        if os.path.isfile(target):
            with open(target) as fh:
                f_targets = fh.readlines()

            for target in f_targets:
                target = target.strip()
                if not target.startswith('#'):
                    targets_list.append(target.upper())

        else:
            targets_list.append(target.upper())

    return targets_list


def parse_models(models):
    """Parser for user input"""
    all_models = cmd.get_names('objects')

    if models.lower() == 'all':
        models = all_models
    else:
        models = models.split(" ")

        for model in models:
            if model not in all_models:
                logging.error("model '{0}' does not exist.".format(model))

    return models


def parse_chains(chains):
    """Parser for user input"""
    all_models = cmd.get_names('objects')
    all_chains = list()

    for model in all_models:
        all_chains.extend(cmd.get_chains(model))

    all_chains = list(set(all_chains))

    if chains.lower() == 'all':
        chains = all_chains

    else:

        chains = [chain.upper() for chain in chains.split(" ")]

        for chain in chains:
            if chain not in all_chains:
                logging.error("chain '{0}' does not exist.".format(chain))

    return chains


def parse_search(search):
    """Parser for user input"""
    if re.match(r'(?:aa|amino|aminoacid)s?', search, re.I):
        search = "aminoacids"
    elif re.match(r'(?:na|nucleic|nucleicacid)s?', search, re.I):
        search = "nucleicacids"
    else:
        logging.error("parameter 'search' is invalid")

    return search


def parse_firstonly(firstonly):
    """Parser for user input"""
    if re.match(r'(?:true|t|1)', firstonly, re.I):
        firstonly = True
    elif re.match(r'(?:false|f|0)', firstonly, re.I):
        firstonly = False
    else:
        logging.error("parameter 'firstonly' is not a valid boolean value")

    return firstonly


def parse_gapcost(gapcost):
    """Parser for user input"""
    try:
        float(gapcost)
    except ValueError:
        logging.error("parameter 'gapcost' is not a valid flaot value")

    return float(gapcost)


def parse_minscore(minscore):
    """Parser for user input"""
    try:
        if float(minscore) > 100 or float(minscore) < 0:
            logging.error("minscore value is not in range of 0 and 100")
    except ValueError:
        logging.error("parameter 'minscore' is not a valid int value")

    return float(minscore)

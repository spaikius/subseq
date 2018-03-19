import re
from pymol import cmd

def get_all_models():
    models = cmd.get_names()

    if not bool(models):
        raise Exception("No models are currently opened")

    return models
  

def get_all_chains(models):
    chain_list = list()
    for model in models:
        for chain in cmd.get_chains(model):
            chain_list.append(chain)

    return list(set(chain_list))


def parse_str_to_list(string):
     # split by any seperator (, . / etc.)
    raw_str = re.sub('([A-Za-z0-9_]+)', r'\1', string)

    # remove symbols '[', ']', ''', '"', white spaces
    raw_str = re.sub('[\[\]\'\"\s]', '', raw_str)

    return raw_str.split(',')

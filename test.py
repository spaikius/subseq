from pymol import cmd
import re
import os


def subseq_re(target=None, chains='all', firstonly=False, models='all'):
	# Check if target is given
	if target is None:
		print('[Error] target is not specified.\n [Help] For more help type: subseq --help')

	# Print manual if target value is [-h, --h, -he, --he, -hel, --hel, -help, --help]
	if re.match(r'^-{1,2}h(?:elp|el|e|)', target, re.I):
		print('help')

	if models.lower() == 'all':
		models = get_all_models()
	else:
		models = parse_str_to_list(models)

	if chains.lower() == 'all':
		chains = get_all_chains()

cmd.extend('subseq.re', subseq_re)
cmd.extend('subseq.la', )


# -----------------------------------------------------------------------------
# HELPER FUNCTIONS

def get_all_models():
    """ return all oppend pymol models """
    
    models = [model for model in cmd.get_names()]

    # if a list is empty -> no models are currently avialable
    if not bool(models):
        raise NoModelsError("No models are currently opened")

    return models
  

def get_all_chains(models):
    """ return a list of chains for given models """

    chain_list = list()

    for model in models:
        for chain in cmd.get_chains(model):
            chain_list.append(chain.upper())

    return list(set(chain_list))


def get_chains(model):
    """ return a list of chains from specific model """

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

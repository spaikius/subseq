"""Description
This module is desgined to extract and access data from pymol.
"""
from pymol import cmd
from HelperFunctions import get_chains


class Data:
	"""
    This class creates an empty dictionary and fills it with data
    from pymol cmd.iterate command.

    self.data schema:
    self.data = {   model_name: {
                        chain_name:{
                            sequence: str
                            ids: [...]
                        },
                        ...
                    },
                    ...
                }
	"""
	# Dictionary to look up the one letter codes
    one_letter = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
        'CA': ''
    }

    def __init__(self, models, chains):
        self.data = None
        self.models = models
        self.chains = chains

        self.construct_empty_data_dict()
        self.fill_empty_data_dict()

    def __getitem__(self, key):
        # 
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                raise Exception("DATA: bad key {}".format(key))
        else:
            return self.data[key]

    def keys(self):
    	"""Returns all self.data top level keys"""
        return self.data.keys()

    def construct_empty_data_dict(self):
    	"""
 		Initializes data structure self.data 
        structure schema:
            model: {chain: {sequence: '', ids: list()}}
            
        model and chain names are in self.models (list) and 
        self.chains (list) respectively
    	"""
        empty_dict = dict()
        sequence = 'sequence'
        ids = 'ids'

        for model in self.models:
            empty_dict[model] = dict()

            for chain in get_chains(model):
                # skip if chain is not requested
                if chain not in self.chains:
                    continue

                empty_dict[model][chain] = dict()
                empty_dict[model][chain][sequence] = ''
                empty_dict[model][chain][ids] = list()

        self.data = empty_dict

    def fill_empty_data_dict(self):
        aa_dict = self.get_data_from_pymol()
        self.replace_to_one_letter(aa_dict)
        self.fill_data(aa_dict)

    @staticmethod
    def get_data_from_pymol():
    	"""
    	 Extracts data from pymol using `cmd.iterate` command.
            Returns a dictionary:
                aa_dict = { 'aa_list': [[resn, resi, chain, model], ...}
                where
                    resn - position number,
                    resi - 3 letter aa code
                    chain - chain name
                    model - model name
    	"""
        aa_dict = dict()
        aa_dict['aa_list'] = list()

        # iterate through all c-alpha atoms and append
        # position number, three letter code, corresponding chain and model
        # to aa_dict['aa_list']
        cmd.iterate("(name ca)",
                    "aa_list.append([resn, resi, chain.upper(), model.upper()])",
                    space=aa_dict)

        return aa_dict

    def replace_to_one_letter(self, aa_dict):
    	"""Replaces all 3 letter aa code to 1 letter aa code"""
        for aa_list in aa_dict['aa_list']:
            aa_list[0] = self.one_letter[aa_list[0]]

    def fill_data(self, aa_dict):
    	"""
        Initializes self.data[model][chain] keys:
    		- sequence: aa chain sequence
    		- ids: list of ids
    	"""
        sequence = 'sequence'
        ids = 'ids'

        for resn, resi, chain, model in aa_dict['aa_list']:
            # Skip if model or chain is not requested
            if model not in self.models or chain not in self.chains:
                continue

            self.data[model][chain][sequence] += resn
            self.data[model][chain][ids].append(resi)

    def filter_data(self):
    	"""Remove all blank attributes in self.data dictionary"""
        for model in self.data.keys():
        	# remove all empty models and warn if found any
            if not self.data[model]:
                print(" Warning: model {} does not contain any of these chains:\n{}"
                      .format(model, self.chains))
                del self.data[model]

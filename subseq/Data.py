"""Description

This module is desgined to extract and access data from pymol.

Data(class):
    Description
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

    Class attributes:
        - one_later(dict): static
            dictionary to look up the one letter codes

        - __init__(models: list, chains: list): typed
            Constructor

        - __getitem__(key: tuple/str) -> dict/list/str: typed
            Getter. Returns self.data inner child

        - keys() -> list: typed
            Returns all self.data top level keys
        
        - construct_empty_data_dict(): typed
            Initializes data structure self.data 
            schema:
                model: {chain: {sequence: '', ids: list()}}
            
            model and chain names are provided in self.models (list) and 
            self.chains (list) respectively

        - get_data_from_pymol() -> dict: static
            Extracts data from pymol using `cmd.iterate` command.
            Returns a dictionary:
                aa_dict = { 'aa_list': [[resn, resi, chain, model], ...}
                where
                    resn - position number,
                    resi - 3 letter aa code
                    chain - chain name
                    model - model name

        - fill_empty_data_dict(aa_dict: dict): typed
            Fills self.data using get_data_from_pymol() constructed dictionary.
            Calls replace_to_one_letter(), fill_data(), filter_data() respectively 
        
        - replace_to_one_letter(aa: dict): typed
            Replaces all 3 letter aa code to 1 letter aa code
            for get_data_from_pymol() constructed dictionary

        - fill_data(): typed
            Initializes self.data[model][chain] keys:
                - sequence: aa chain sequence
                - ids: list of ids
            please note: 

        - filter_data(): typed
            Remove all blank attributes in self.data dictionary


    Accessing data through class object:
        class_obj['model_name', 'chain_name', 'sequence']
        or
        class_obg['model_name']['chain_name']['sequence']

"""
from pymol import cmd
from HelperFunctions import get_chains


class Data:
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
        return self.data.keys()

    def construct_empty_data_dict(self):
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
        aa_dict = dict()
        aa_dict['aa_list'] = list()

        # iterate through all c-alpha atoms and append
        # position number, three letter code, corresponding chain and model
        # to aa_dict['aa_list']
        cmd.iterate("(name ca)",
                    "aa_list.append([resn, resi, chain, model])",
                    space=aa_dict)

        return aa_dict

    def replace_to_one_letter(self, aa_dict):
        for aa_list in aa_dict['aa_list']:
            aa_list[0] = self.one_letter[aa_list[0]]

    def fill_data(self, aa_dict):
        sequence = 'sequence'
        ids = 'ids'

        for resn, resi, chain, model in aa_dict['aa_list']:
            # Skip if model or chain is not requested
            if model not in self.models or chain not in self.chains:
                continue

            self.data[model][chain][sequence] += resn
            self.data[model][chain][ids].append(resi)

    def filter_data(self):
        # remove all empty models and warn if found any
        for model in self.data.keys():
            if not self.data[model]:
                print(" Warning: model {} does not contain any of these chains:\n{}"
                      .format(model, self.chains))
                del self.data[model]

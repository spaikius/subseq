"""Description
This module is desgined to extract and access data from pymol.
"""
from pymol import cmd


class Data:
    """
    This class is designed to extract data from pymol using cmd.iterate command
    and makes it accessible through class object.

    self.data schema:
    self.data = {
                    model: {
                        chain:{
                            sequence: str
                            ids: list
                        },
                        ...
                    },
                    ...
                }
    """

    # Dictionary to look up the one letter codes
    aa_one_letter = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    }

    na_one_letter = {
        'DG': 'G', 'DA': 'A', 'DT': 'T', 'DC': 'C', 'DU': 'U',
        'G': 'G', 'A': 'A', 'T': 'T', 'C': 'C', 'U': 'U'
    }

    def __init__(self, models, chains, search_for, replace_with):
        self.data = None
        self.models = models
        self.chains = chains
        self.search_for = search_for.lower()
        self.replace_with = replace_with

        self.construct_data_dict()
        self.fill_data_dict()

    def __getitem__(self, key):
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                raise KeyError("bad key {}".format(key))
        else:
            return self.data[key]

    def keys(self):
        """Returns all self.data top level keys"""
        return self.data.keys()

    def construct_data_dict(self):
        """
        Initializes data structure self.data_dict
        structure schema:
            model: {chain: {sequence: string, ids: list}}
        """
        data_dict = dict()

        for model in self.models:
            data_dict[model] = dict()

            for chain in cmd.get_chains(model):
                # skip if chain is not requested
                if chain not in self.chains:
                    continue

                data_dict[model][chain] = dict()
                data_dict[model][chain]['sequence'] = ''
                data_dict[model][chain]['ids'] = list()

        self.data = data_dict

    def fill_data_dict(self):
        atoms_dict = self.get_data_from_pymol()

        if self.search_for == 'aminoacids':
            self.replace_to_aa_one_letter(atoms_dict)
        elif self.search_for == 'nucleicacids':
            self.replace_to_na_one_letter(atoms_dict)

        self.fill_data(atoms_dict)
        self.filter_data()

    def get_data_from_pymol(self):
        """
         Extracts data from pymol using `cmd.iterate` command.
            Returns a dictionary:
                atoms_dict = { 'main_atoms': [[resn, resi, chain, model], ...}
                where
                    resn - position number,
                    resi - 3 letter aa code
                    chain - chain name
                    model - model name
        """
        atoms_dict = dict()
        atoms_dict['main_atoms'] = list()
        # iterate through all c-alpha atoms and append
        # position number, three letter code, corresponding chain and model
        # to atoms_dict['main_atoms']
        if self.search_for == 'aminoacids':
            cmd.iterate("(name ca)",
                        "main_atoms.append([resn, resi, chain, model])",
                        space=atoms_dict)

        if self.search_for == 'nucleicacids':
            cmd.iterate("(resn G+C+A+T+U+DG+DC+DA+DT+DU and name C1')",
                        "main_atoms.append([resn, resi, chain, model])",
                        space=atoms_dict)

        return atoms_dict

    def replace_to_aa_one_letter(self, atoms_dict):
        """Replaces all 3 letter amino acid code to 1 letter aa code"""
        for main_atoms in atoms_dict['main_atoms']:
            try:
                main_atoms[0] = self.aa_one_letter[main_atoms[0]]
            except KeyError:
                main_atoms[0] = self.replace_with

    def replace_to_na_one_letter(self, atoms_dict):
        """Replaces all DNA or RNA nucleic acid code to 1 letter aa code"""
        for main_atoms in atoms_dict['main_atoms']:
            try:
                main_atoms[0] = self.na_one_letter[main_atoms[0]]
            except KeyError:
                main_atoms[0] = self.replace_with

    def fill_data(self, atoms_dict):
        """
        Initializes self.data[model][chain] keys:
            - sequence: aa chain sequence
            - ids: list of ids
        """
        for resn, resi, chain, model in atoms_dict['main_atoms']:
            # Skip if model or chain is not requested
            if model not in self.models or chain not in self.chains:
                continue

            self.data[model][chain]['sequence'] += resn
            self.data[model][chain]['ids'].append(resi)

    def filter_data(self):
        for model in self.data.keys():
            for chain in self.data[model].keys():
                if self.data[model][chain]['sequence'] == '':
                    self.data[model].pop(chain)

        for model in self.data.keys():
            if len(self.data[model].keys()) == 0:
                self.data.pop(model)
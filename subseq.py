from pymol import cmd
import re
import os
from random import randint


class InvalidPairError(Exception):
    """ Invalid key pair Exception """
    pass

class NoModelsError(Exception):
    """ No models found Exception """
    pass

class BadParameterError(Exception):
    """ Bad Parameter Exception """
    pass

class BadRegExSyntaxError(Exception):
    """ Regular Expresion syntax Exception """
    pass

class InvalidMatrixFormatError(Exception):
    """ Invalid substitution matrix Exception """
    pass

# MAIN
def subseq(*argv, **kwargs):
    # If user just typed 'subseq'
    # kwargs contains `_self` key, so it's length is 1, not 0
    if len(argv) == 0 and len(kwargs) == 1:
        print("For help type: subseq help")
        return

    # Print usage manual if user typed `subseq help` 
    if 'help' in argv:
       print(usage_message)
       return

    # Print warrning message if user forgot to seperate kwargs with comma.
    # Or targets value is not within parentheses if RegExp quantifier `{n,m}` is used
    if len(kwargs) != 1 and len(argv) > 0:
        print(bad_arguments_message)
        return

    try:
        parameters = ParseKwargs(kwargs)

        data = Data(parameters['models'], parameters['chains'])

        search_result = None
        if parameters['method'] == 're':
            search_result = subseq_re(parameters['target'], data)
        elif parameters['method'] == 'la':
            search_result = subseq_la(  parameters['target']
                                      , data
                                      , parameters['submatrix']
                                      , parameters['gapcost']
                                      , parameters['minscore'])

        if search_result is not None:
            select(search_result)
        else:
            print("Nothing can be found for given target: {0}".format(parameters['target']))

    except (  BadParameterError
            , InvalidPairError
            , BadRegExSyntaxError
            , NoModelsError
            , InvalidMatrixFormatError
            , InvalidPairError ) as e:
        print(' [Error] ' + str(e))

    except Exception:
        # Unknown exception occured.
        # Reraise exception and let Pymol to handle it. 
        raise

cmd.extend('subseq', subseq)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS

def get_all_models():
    """ return all oppend pymol models """
    
    models = [model.upper() for model in cmd.get_names()]

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

# END OF HELPER FUNCTIONS
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CLASS PARSEKWARGS
class ParseKwargs:
    """
    This class is for parsing and validating raw kwargs (keyword arguments)
    from Pymol command line.
    """
    method_types = [
        're',  # Regular Expression (default)
        'la'   # Local alignment
    ]

    default_method = method_types[0]
    default_gap_cost = 10
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
        if self.parameters['method'] is 're':
            return

        if self.parameters['submatrix'] is None:
            raise BadParameterError("No substitution matrix file")
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

# END OF CLASS PARSEKWARGS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CLASS DATA
class Data:
    """
    This class is desgined to extract data from pymol cmd.iterate command
    and makes it accessible through class object.

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
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                raise InvalidPairError("bad key {}".format(key))
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

# END OF CLASS DATA
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CLASS SUBMATRIX

class SubMatrix:
    def __init__(self, matrix_path):
        self.matrix = None
        self.name = os.path.basename(matrix_path)
        self.load_matrix(matrix_path)

    def load_matrix(self, matrix_path):
        with open(matrix_path, 'r') as fh:
            matrix = fh.read()

        lines = matrix.strip().split('\n')
        # remove comments
        lines = [line for line in lines if line[0] != '#']

        header = lines.pop(0)
        columns = header.split()
        matrix = dict()

        for row in lines:
            entries = row.split()
            row_name = entries.pop(0)
            matrix[row_name] = dict()

            if len(entries) != len(columns):
                raise InvalidMatrixFormatError('columns and rows counts does not match\n file: {}'
                                .format(self.matrix))
            for column_name in columns:
                matrix[row_name][column_name] = entries.pop(0)

        self.matrix = matrix

    def __getitem__(self, key): 
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                 raise InvalidPairError("Bad key pair: {}".format(key))
        else:
            return self.matrix[key]

    def get_name(self):
        return self.name

# END OF CLASS SUBMATRIX
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LOCAL ALIGNMENT

def subseq_la(target, data, matrix, gap_cost, min_score):
    """
    
    """
    sub_matrix = SubMatrix(matrix)
    max_score = calculate_max_score(target, sub_matrix)
    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            SW = SmithWaterman(target, sequence, gap_cost, sub_matrix)

            # Skip if alignment best score is less than minimum passing score
            if (float(SW.get_best_score()) / max_score) * 100 < min_score:
                continue

            for i, j in SW.get_coordinates():
                aligned_target, aligned_sequence, start_i, start_j = SW.get_traceback(i, j)

                # start_j - 1, because start_j corresponds to index from Smith-Waterman score table
                start_id = data[model][chain]['ids'][start_j - 1]

                # end = start index + alignment length (without gaps)
                subject_end = len([i for i in aligned_sequence if i != '-']) + start_j - 1 
                
                # subject_end - 1, because list index starts at 0               
                end_id = data[model][chain]['ids'][subject_end - 1]

                match_list.append((model, chain, start_id, end_id))

                alignment_string, identities, gaps, mismatches = \
                    create_alignment_string(aligned_target, aligned_sequence)

                print_alignment(  model, chain, target, sequence, sub_matrix.get_name(), gap_cost
                                , SW.get_best_score(), max_score, identities, mismatches, gaps
                                , aligned_target, aligned_sequence, alignment_string
                                , start_i, start_j, data[model][chain]['ids'])


    return match_list if len(match_list) != 0 else None


def calculate_max_score(target, sub_matrix):
    """ calculate the best possible alignemnt score for given target """
    max_score = 0
    for aa in target:
        max_score += int(sub_matrix[aa, aa])

    return max_score


def create_alignment_string(aligned_seq1, aligned_seq2):
    """
    Constructs alignment string for both aligned sequences

    KTGTA
    :| :|  <-- alignment string
    PT-KA

    where ':' - mismatch, ' ' - gap, '|' - match
    
    returns alignment string, match count, mismatch count, gap count
    """
    identities, gaps, mismatches = 0, 0, 0
    alignment_string = ''

    for aa1, aa2 in zip(aligned_seq1, aligned_seq2):
        if aa1 == aa2:
            alignment_string += '|'
            identities += 1

        elif '-' in (aa1, aa2):
            alignment_string += ' '
            gaps += 1

        else:
            alignment_string += ':'
            mismatches += 1

    return alignment_string, identities, gaps, mismatches


def print_alignment(  model, chain, target, sequence, substitution_matrix_name, gap_cost
                    , alignment_score, max_score, identities, mismatches, gaps
                    , aligned_target, aligned_sequence, alignment_string
                    , target_start, subject_start, ids_list):

    """
    Prints BLAST like alginment for both sequences
    """

    a_len = len(aligned_sequence)

    print("\n")
    print("Model: {0}, chain: {1}".format(model, chain))
    print("Target length: {0} {1}".format(len(target), target[0:40]))
    print("Subject length: {0} {1}".format(len(sequence), sequence[0:40]))
    print("Substitution matrix: {0}".format(substitution_matrix_name))
    print("Gap cost: {0}".format(gap_cost))
    print("\n")
    print("Alignment score: {0}/{1} ({2:.1%})"
        .format(alignment_score, max_score, float(alignment_score) / max_score))

    print("Identities: {0}/{1} ({2:.1%})"
        .format(identities, a_len, float(identities) / a_len))

    print("Mismatches: {0}/{1} ({2:.1%})"
        .format(mismatches, a_len, float(mismatches) / a_len))

    print("Gaps:       {0}/{1} ({2:.1%})"
        .format(gaps, a_len, float(gaps) / a_len))

    for i in range(0, a_len, 60):
        target_slice = aligned_target[i: i + 60]
        subject_slice = aligned_sequence[i: i + 60]
        alignment_slice = alignment_string[i: i + 60]

        target_start = i + target_start
        subject_start = i + int(ids_list[subject_start - 1])

        target_end = len([i for i in target_slice if i != '-']) + target_start - 1
        subject_end = len([i for i in subject_slice if i != '-']) + subject_start - 1

        print("Target  {0:<4} {1} {2}".format(target_start, target_slice, target_end))
        print(' ' * 13 + "{0}".format(alignment_slice))
        print("Subject {0:<4} {1} {2}".format(subject_start, subject_slice, subject_end))
        print("\n")

    print('-' * 60)

# END OF LOCAL ALIGNMENT
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CLASS SMITHWATERMAN
class SmithWaterman:
    """
    This class performs nucleotide or protein sequence (depending on given 
    substitution matrix) alignment using the Smith-waterman algorithm
    """
    def __init__(self, target, sequence, gap_cost, sub_matrix):
        self.target = [i for i in target]
        self.sequence = [i for i in sequence]
        self.gap_cost = gap_cost
        self.sub_matrix = sub_matrix

        self.best_score = 0
        self.best_score_coordinates = list()
        
        """
        Initialize score matrix
                 S  E  Q  U  E  N  C  E
            [[0, 0, 0, 0, 0, 0, 0, 0, 0],
          T  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          A  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          R  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          G  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          E  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          T  [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        """
        self.score_matrix = [[0 for i in range(len(sequence) + 1)] for j in range(len(target) + 1)] 

        self.fill_score_matrix()

    def __getitem__(self, i, j): 
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                 raise InvalidPairError("Bad key pair: {}".format(key))
        else:
            return self.matrix[key]

    def get_coordinates(self):
        return self.best_score_coordinates

    def get_best_score(self):
        return self.best_score

    def get_traceback(self, i, j):
        """
        Finds the optimal path through the self.score_matrix.
        Returns constructed alignment strings and coordinates where alignment begins
        """
        aligned_target = list()
        aligned_subject = list()

        end, diag, up, left = range(4)

        move = self.next_move(i, j)

        while move != end:
            if move == diag:
                aligned_target.append(self.target[i - 1])
                aligned_subject.append(self.sequence[j - 1])

                i -= 1
                j -= 1

            elif move == up:
                aligned_target.append(self.target[i - 1])
                aligned_subject.append('-')

                i -= 1

            elif move == left:
                aligned_target.append('-')
                aligned_subject.append(self.sequence[j - 1])

                j -= 1

            move = self.next_move(i, j)

        aligned_target.append(self.target[i - 1])
        aligned_subject.append(self.sequence[j - 1])

        aligned_target = ''.join(reversed(aligned_target))
        aligned_subject = ''.join(reversed(aligned_subject))

        return aligned_target, aligned_subject, i, j

    def next_move(self, i, j):
        """
        Looks for the next move during traceback.
        Moves are determined by the score of three upper-left, left and up 
        in the self.score_matrix elements
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]
        achieved_score = self.score_matrix[i][j]
        diag = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if achieved_score == diag + int(self.sub_matrix[aa1, aa2]):
            # return diagnol move if diagnol move is greater than 0 else return END
            return 1 if diag > 0 else 0

        if achieved_score == up - self.gap_cost:
            # return up move if up move is greater than 0 else return END
            return 2 if up > 0 else 0

        if achieved_score == left - self.gap_cost:
            # return left move if left move is greater than 0 else return END
            return 3 if left > 0 else 0

    def fill_score_matrix(self):
        """
        Fills self.score_matrix with scores representing trial 
        alignments of the two sequences
        """
        for i in range(1, len(self.score_matrix)):
            for j in range(1, len(self.score_matrix[i])):
                score = self.calculate_score(i, j)

                self.score_matrix[i][j] = score

                if score > self.best_score:
                    self.best_score = score
                    self.best_score_coordinates = [(i, j)]
                elif score == self.best_score:
                    self.best_score_coordinates.append((i, j))

    def calculate_score(self, i, j):
        """
        Calculates score for given i and j position in the self.score_matrix.
        The score is based on the upper-left, left and up elements in self.score_matrix
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]

        similarity = self.sub_matrix[aa1, aa2]

        diag_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(0, diag_score, up_score, left_score)

# END OF CLASS SMITHWATERMAN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# REGULAR EXPRESSION SEARCH

def subseq_re(target, data):
    """
    workflow:
        1) create a RegExp object from target (function parameter)
        2) scan data (function parameter) by using RegExp object
        3) append information about match to match_list
        4) return match_list if its length is not 0 else return None
    """
    sequence = 'sequence'
    ids = 'ids'
    match_list = list()

    # RegExp validation
    try:
        # re.I - ignore case sensetive
        re_target = re.compile(target, re.I)
    except:
        raise BadRegExSyntaxError('Bad syntax - target(RegExp): ' + str(target))

    # scan data by using RegExp object
    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain][sequence]):

                start_id = data[model][chain][ids][match.start()]
                end_id = data[model][chain][ids][match.end() - 1]
                
                match_list.append((model, chain, start_id, end_id))

    # return match list if its length is not 0 else return None
    return match_list if len(match_list) != 0 else None

# END OF REGULAR EXPRESSION SEARCH
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# SELECT

def select(select_list):
    lower_bound = 0
    upper_bound = 100000
    # random select ID
    select_id = 'ss_' + str(randint(lower_bound, upper_bound))
    # empty select
    cmd.select(select_id, None)

    select_query = None

    for sele_tuple in select_list:
        model = sele_tuple[0]
        chain = sele_tuple[1]
        start_id = sele_tuple[2]
        end_id = sele_tuple[3]

        # select /model/?/chain/resi-resi
        select_query = " | /{0}//{1}/{2}-{3}".format(model, chain, start_id, end_id)

        # Execute and append selection to select_id
        cmd.select(select_id, select_id + select_query)

# END OF SELECT
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PRINT MESSAGES
usage_message = """
Usage: subseq target=<str>, method=<str>, submatrix=<str>, models=<list>
              , chains=<list>, gapcost=<float>, minscore=<float>

Exaple usage: subseq target=KTGT, method=la, chains=[A, B, T], submatrix=PATH/TO/MATRIX

!!! Important !!!
Please note: each keyword parameter should be seperated with comma (,)
Please note: target value should be within parentheses if quantifier {n,m} is used in Regular Expressions.
             Example: target=(GT{3,})
                             ^      ^


Parameters:
    help                            ; Prints usage manual

Keyword parameters:
    target=<str>        Required    ; Target sequence
                                      Examples:
                                       If method type is re:
                                        - target=KTGTAVU
                                        - target=(TATA.{3,5}ATG(.{3,4}){3,})
                                       If method type is la:
                                        - target=KTGAT


    method=<str>        Optional    ; Method search type.
                                      - 're' for Regular Expression
                                      - 'la' for local alignment. Smith-Waterman 
                                      Default value: 're'

    models=<list>       Optional    ; The list of models that will be used for target search.
                                      If the list is not provided, then search for target will be 
                                      performed in all available models
                                      Example: models=[5ara, 2cif, a4s2]

    chains=<list>       Optional    ; The list of chains that will be used for target search.
                                      If the list is not priveded, then search for target will be
                                      performed in all available model chains
                                      Example: chains=[A, T, X, Q]

    submatrix=<PATH>    Required    ; Path to substitution matrix for local alignment

    gapcost=<float>     Optional    ; The linear gap cost for local alignment
                                      Default value: 10

    minscore=<float>    Optional    ; The minimum score in precentages for throwing off low score alignments
                                      in local alignment search.
                                      Default value: 51.00 ( 51% )
                                      Example: minscore=75.25
"""

bad_arguments_message = """
Found nonseperated keyword arguments or target= value is not within parentheses.

Please check if all keyword arguments are seperated with comma (,) or target= value
is within parentheses if Regular Expression quantifier is used.

For help type: subseq help
"""
# END OF PRINT MESSAGES
# -----------------------------------------------------------------------------

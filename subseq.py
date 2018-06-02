import re
import os
import logging
import string

from pymol import cmd, stored

__author__ = "Rimvydas Noreika"
__credits__ = ["Justas Dapkūnas"]
__version__ = "1.0.1"
__status__ = "Development"


def subseq_re(
        targets='', chains='all', search='aminoAcids', firstonly='False',
        models='all', sele='ss-{method}-{id}-{target}'):
    """
DESCRIPTION
    subseq - tool for searching target sequences using Regular Expressions

USAGE
    subseq targets=<list|FILE>, chains=<list>, search=<str>, 
           firstonly=<bool>,    models=<list>, sele=<str>

    Example usage: subseq KTGT, [A, B, C], firstonly=True, search=nucleicAcids

IMPORTANT
    All modified amino or nucleic acids are replaced with: X

    If quantifier {n,m} is used in Regular Expressions then target value should
    be within parentheses, single or double quatation marks
             
    Example: target= (GT{3,})   or  'GT{3,}'   or  "GT{3,}"
                     ^      ^       ^      ^       ^      ^
PARAMETERS
    targets=<list|FILE>     ; Target sequence
                              Examples:
                                - targets=KTGTAVU
                                - targets="TATA.{3,5}ATG(.{3,4}){3,}"
                                - targets=[SIS, KATK, "AK{3,4}"]
                                - targets=PATH/TO/TARGETS_FILE
    
    chains=<list>           ; The list of chains 
                              Examples: 
                                - chains=A
                                - chains=[A, AT, X, Q]
                              Default: all
    
    search=<str>            ; Search for nucleic acids or amino acids sequence
                                - for amino acids: aminoAcids, amino, aa
                                - for nucleic acids: nucleicAcids, nucleic, na
                              Default value: aminoAcids

    firstonly=<bool>        ; If firstonly is False (0) then select all matches
                              If firstonly is True  (1) then select first match
                              Default: False

    models=<list>           ; The list of models
                              Examples:
                                - models=5ara
                                - models=[5ara, 2cif, a4s2]
                              Default: all

    sele=<str>              ; Selection name
                              Tokens:
                                - {method} - used method for sequence search
                                - {target} - first 10 targets alphabet letters
                                - {id}     - id
                              Default: 'ss-{method}-{id}-{target}'

SEE ALSO
    subseq.local, subseq.global

SUBSEQ                          2018-06-01
    """
    
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.error = CallCounter(logging.error)

    targets = parse_targets(targets)
    chains = parse_chains(chains)
    search = parse_search(search)
    firstonly = parse_firstonly(firstonly)
    models = parse_models(models)

    if logging.error.counter is not 0:
        logging.info("{0} errors were found. ".format(logging.error.counter) 
            + "Please see above messages for more information")

        return
        

    data = Data(models, chains, search, replace_with='X')
    
    search_results = None

    for target in targets:
        try:
            search_results = subseq_re_search(target, data, firstonly, search)

        except Exception as e:
            logging.warning("RegExp for {0}: {1}".format(target, e))
            continue

        if search_results is not None:
            select(search_results, target, sele, method='re')

        else:
            logging.info("Nothing can be found for given target: {0}"
                .format(target))


def subseq_local_alignment(
        targets='', submatrix='blossum62', chains='all', search='aminoacids', 
        firstonly='False', gapcost='10', minscore='51', models='all',
        sele='ss-{method}-{id}-{target}'):
    """
DESCRIPTION
    subseq.local - tool for searching target sequences using local alignment

USAGE
    subseq.local targets=<list|FILE>, submatrix=<FILE> chains=<list>,
                 search=<str>, firstonly=<bool>, gapcost=<float>,
                 minscore=<float>, models=<list>, sele=<str>

    Example usage: subseq.local KTGT, blossum62, firstonly=True, gapcost=12.5

IMPORTANT
    All modified amino or nucleic acids are replaced with: X

PARAMETERS
    targets=<list|FILE>     ; Target sequence
                              Examples:
                                - targets=KTGTAVU
                                - targets="TATA.{3,5}ATG(.{3,4}){3,}"
                                - targets=[SIS, KATK, "AK{3,4}"]
                                - targets=PATH/TO/TARGETS_FILE
    
    submatrix=<FILE>        ; Path to substitution matrix file
                              Default: blossum62

    chains=<list>           ; The list of chains 
                              Examples: 
                                - chains=A
                                - chains=[A, AT, X, Q]
                              Default: all
    
    search=<str>            ; Search for nucleic acids or amino acids sequence
                                - for amino acids: aminoAcids, amino, aa
                                - for nucleic acids: nucleicAcids, nucleic, na
                              Default value: aminoAcids

    firstonly=<bool>        ; If firstonly is False (0) then select all matches
                              If firstonly is True  (1) then select first match
                              Default: False

    gapcost=<float>         ; The linear gap cost for local alignment
                              Default value: 10

    minscore=<float>        ; The minimum score in precentages for throwing off
                              low score alignments in alignment search
                              Example: minscore=75.25
                              Default value: 51.00 ( 51% ) 

    models=<list>           ; The list of models
                              Examples:
                                - models=5ara
                                - models=[5ara, 2cif, a4s2]
                              Default: all

    sele=<str>              ; Selection name
                              Tokens:
                                - {method} - used method for sequence search
                                - {target} - first 10 targets alphabet letters
                                - {id}     - id
                              Default: 'ss-{method}-{id}-{target}'

SEE ALSO
    subseq, subseq.global

SUBSEQ                          2018-06-01
    """
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.error = CallCounter(logging.error)

    targets = parse_targets(targets)
    chains = parse_chains(chains)
    search = parse_search(search)
    firstonly = parse_firstonly(firstonly)
    gapcost = parse_gapcost(gapcost)
    minscore = parse_minscore(minscore)
    models = parse_models(models)

    if logging.error.counter is not 0:
        logging.info("{0} errors were found. ".format(logging.error.counter) 
            + "Please see above messages for more information")

        return
    
    if search is 'nucleicacids' and submatrix is 'blossum62':
        submatrix = 'nucleicmatrix'

    data = Data(models, chains, search, replace_with='X')

    for target in targets:
        try:
            search_results = subseq_la_search(target, data, submatrix, 
                                          gapcost, minscore, firstonly)
        except Exception as e:
            logging.error("{0}".format(e))
            continue

        if search_results is not None:
            select(search_results, target, sele, method='local')

        else:
            logging.info("Nothing can be found for given target: {0}"
                .format(target))


def subseq_global_alignment():
    """UNDER DEVELOPMENT"""
    pass

cmd.extend('subseq', subseq_re)
cmd.extend('subseq.local', subseq_local_alignment)
cmd.extend('subseq.global', subseq_global_alignment)
stored.id = 0

class CallCounter:
    """Decorator to determine number of calls for a method"""

    def __init__(self, method):
        self.method = method
        self.counter = 0

    def __call__(self, *args, **kwargs):
        self.counter += 1
        return self.method(*args, **kwargs)


def parse_targets(targets):
   
    if targets is '':
        logging.error("parameter 'targets' is not specified.")
        return

    if os.path.isfile(targets):
        with open(targets) as targets_file:
            targets = targets_file.readlines()

        targets = [ target.strip() for target in targets
                        # Skip comments 
                        if not target.startswith('#') 
                  ] 
    else:
        targets = parse_str_to_list(targets)

    return targets


def parse_models(models):
    all_models = cmd.get_names('objects')

    if models.lower() == 'all':
        models = all_models
    else:
        models = parse_str_to_list(models)

        for model in models:
            if model not in all_models:
                logging.error("model '{0}' does not exist.".format(model))

    return models


def parse_chains(chains):
    
    all_models = cmd.get_names('objects')
    all_chains = list()

    for model in all_models:
        all_chains.extend(cmd.get_chains(model))

    all_chains = list(set(all_chains))

    if chains.lower() == 'all':
        chains = all_chains

    else:        

        chains = [ chain.upper() for chain in parse_str_to_list(chains) ]

        for chain in chains:
            if chain not in all_chains:
                logging.error("chain '{0}' does not exist.".format(chain))

    return chains


def parse_search(search):
    
    if re.match(r'(?:aa|amino|aminoacid)s?', search, re.I):
        search = "aminoacids"
    elif re.match(r'(?:na|nucleic|nucleicacid)s?', search, re.I):
        search = "nucleicacids"
    else:
        logging.error("parameter 'search' is invalid")

    return search


def parse_firstonly(firstonly):

    if re.match(r'(?:true|t|1)', firstonly, re.I):
        firstonly = True
    elif re.match(r'(?:false|f|0)', firstonly, re.I):
        firstonly = False
    else:
        logging.error("parameter 'firstonly' is not a valid boolean value")        

    return firstonly


def parse_gapcost(gapcost):
    try:
        float(gapcost)
    except ValueError:
        logging.error("parameter 'gapcost' is not a valid flaot value")

    return float(gapcost)


def parse_minscore(minscore):
    try:
        if float(minscore) > 100 or float(minscore) < 0:
            logging.error("minscore value is not in range of 0 and 100")
    except ValueError:
        logging.error("parameter 'minscore' is not a valid int value")

    return float(minscore)


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
        self.fill_empty_data_dict()

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

    def fill_empty_data_dict(self):
        atoms_dict = self.get_data_from_pymol()

        if self.search_for == 'aminoacids':
            self.replace_to_aa_one_letter(atoms_dict)
        elif self.search_for == 'nucleicacids':
            self.replace_to_na_one_letter(atoms_dict)

        self.fill_data(atoms_dict)
        # self.filter_data()

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


def subseq_re_search(target, data, first_only, search_for):
    """
    work flow:
        1) create a RegExp object
        2) scan data by using RegExp object
        3) append information about match to match_list
        4) return match_list if its length is not 0 else return None
    """
    match_list = list()

    target = target.strip("'()\"")

    # Replace wildcards
    if search_for == 'nucleicacids':
        target = target.replace('R', '[AG]')
        target = target.replace('Y', '[CT]')
        target = target.replace('S', '[GC]')
        target = target.replace('W', '[AT]')
        target = target.replace('K', '[GT]')
        target = target.replace('M', '[AC]')
        target = target.replace('B', '[CGT]')
        target = target.replace('D', '[AGT]')
        target = target.replace('H', '[ACT]')
        target = target.replace('N', '.')

    # RegExp validation
    try:
        # re.I - ignore case sensitive
        re_target = re.compile(target, re.I)
    except Exception:
        raise

    # scan data by using RegExp object
    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain]['sequence']):
                start_pos = match.start()

                for _ in range(0, len(match.group())):
                    resi = data[model][chain]['ids'][start_pos]

                    match_list.append((model, chain, resi))

                    start_pos += 1

                if first_only:
                    break
            else:
                continue
            break

    return match_list if len(match_list) else None


def subseq_la_search(target, data, matrix, gap_cost, min_score, first_only):
    """
    """
    # Substitution matrix
    sub_matrix = SubMatrix(matrix)
    
    # The maximum score for given target
    max_score = calculate_max_score(target, sub_matrix)
    
    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            sw = SmithWaterman(target, sequence, gap_cost, sub_matrix)

            # Skip if alignment best score is less than minimum passing score
            if float(sw.get_best_score()) / max_score * 100 < min_score:
                continue

            for i, j in sw.get_coordinates():

                aligned_target, aligned_sequence, start_i, start_j =\
                    sw.get_traceback(i, j)

                start_pos = start_j - 1
                
                for _ in range(0, len(aligned_sequence.replace('-', ''))):
                    resi = data[model][chain]['ids'][start_pos]
                    match_list.append((model, chain, resi))
                    start_pos += 1

                alignment_string, identities, gaps, mismatches = \
                    create_alignment_string(aligned_target, aligned_sequence)

                print_alignment(
                    model, chain, target, sequence, sub_matrix.get_name(),
                    gap_cost, sw.get_best_score(), max_score, identities,
                    mismatches, gaps, aligned_target, aligned_sequence,
                    alignment_string, start_i, start_j,
                    data[model][chain]['ids'])

                if first_only:
                    break
            else:
                continue
            break

    return match_list if len(match_list) != 0 else None


def calculate_max_score(target, sub_matrix):
    """ calculate the best possible alignment score for given target """
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


def print_alignment(
        model, chain, target, sequence, substitution_matrix_name, gap_cost,
        alignment_score, max_score, identities, mismatches, gaps,
        aligned_target, aligned_sequence, alignment_string,
        target_start, subject_start, ids_list):
    """Prints BLAST like alignment"""

    a_len = len(aligned_sequence)

    print("\n")
    print("Model: {0}, chain: {1}".format(model, chain))
    print("Target length: {0} {1}".format(len(target), target[0:40]))
    print("Subject length: {0} {1}".format(len(sequence), sequence[0:40]))
    print("Substitution matrix: {0}".format(substitution_matrix_name))
    print("Gap cost: {0}".format(gap_cost))
    print("\n")
    print("Alignment score: {0}/{1} ({2:.1%})"
          .format(alignment_score, max_score, 
            float(alignment_score) / max_score))

    print("Identities: {0}/{1} ({2:.1%})"
          .format(identities, a_len, float(identities) / a_len))

    print("Mismatches: {0}/{1} ({2:.1%})"
          .format(mismatches, a_len, float(mismatches) / a_len))

    print("Gaps:       {0}/{1} ({2:.1%})"
          .format(gaps, a_len, gaps / a_len))

    for i in range(0, a_len, 60):
        target_slice = aligned_target[i: i + 60]
        subject_slice = aligned_sequence[i: i + 60]
        alignment_slice = alignment_string[i: i + 60]

        target_start = i + target_start
        subject_start = i + int(ids_list[subject_start - 1])

        target_end = len([i for i in target_slice if i != '-']) \
                     + target_start - 1
        
        subject_end = len([i for i in subject_slice if i != '-']) \
                     + subject_start - 1

        print("Target  {0:<4} {1} {2}"
            .format(target_start, target_slice, target_end))
        
        print(' ' * 13 + "{0}".format(alignment_slice))
        
        print("Subject {0:<4} {1} {2}"
            .format(subject_start, subject_slice, subject_end))
        
        print("\n")

    print('-' * 60)


class SmithWaterman:
    """
    This class performs nucleotide or protein sequence (depending on given
    substitution matrix) alignment using the Smith-waterman algorithm
    """

    def __init__(self, target, sequence, gap_cost, sub_matrix):
        self.target = [i for i in target]
        self.sequence = [i for i in sequence]
        self.gap_cost = float(gap_cost)
        self.sub_matrix = sub_matrix

        self.best_score = 0
        self.best_score_coordinates = list()

        '''
        Initialize score matrix
                 S  E  Q  U  E  N  C  E
            [[0, 0, 0, 0, 0, 0, 0, 0, 0],
          T  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          A  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          R  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          G  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          E  [0, 0, 0, 0, 0, 0, 0, 0, 0],
          T  [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        '''
        self.score_matrix = [[0 for _ in range(len(sequence) + 1)] \
                                for _ in range(len(target) + 1)]

        self.fill_score_matrix()

    def __getitem__(self, key):
        if isinstance(key, tuple):
            try:
                x = self
                for k in key:
                    x = x[k]
                return x
            except:
                raise KeyError("Bad key pair: {}".format(key))
        else:
            return self.score_matrix[key]

    def get_coordinates(self):
        """Retruns a list of tuples (i, j) 
        where i and j are coordinates of the best score
        """
        return self.best_score_coordinates

    def get_best_score(self):
        """Returns the best score"""
        return self.best_score

    def get_traceback(self, i, j):
        """Finds the optimal path through the score matrix.
        Returns constructed alignment strings for target and subject and
        values of i, j where alignment begins
        """
        aligned_target = list()
        aligned_subject = list()

        end, diagonal, up, left = range(4)

        move = self.next_move(i, j)

        while move != end:
            if move == diagonal:
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
        """Looks for the next move during traceback.
        Moves are determined by the score of three upper-left, left and up
        in the score matrix
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]
        achieved_score = self.score_matrix[i][j]
        diagonal = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if achieved_score == diagonal + int(self.sub_matrix[aa1, aa2]):
            # return diagonal move or END
            return 1 if diagonal > 0 else 0

        if achieved_score == up - self.gap_cost:
            # return up move or END
            return 2 if up > 0 else 0

        if achieved_score == left - self.gap_cost:
            # return left move or END
            return 3 if left > 0 else 0

    def fill_score_matrix(self):
        """Fills self.score_matrix with scores representing trial alignments
        of the two sequences
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
        """Calculates score for given i and j position in the score matrix
        The score is based on the upper-left, left and up elements
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]

        similarity = self.sub_matrix[aa1, aa2]

        diagonal_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(0, diagonal_score, up_score, left_score)


class SubMatrix:
    """

    """

    def __init__(self, matrix_path):
        self.matrix = None
        self.name = matrix_path
        self.load_matrix(matrix_path)

    def load_matrix(self, matrix_path):
        """Loads substitution matrix"""
        try:
            with open(matrix_path, 'r') as fh:
                matrix = fh.read()
        except Exception:
            matrix = self.matrices(matrix_path)

            if matrix is None:
                raise

        lines = matrix.strip().split('\n')
        # remove comments
        lines = [line for line in lines if not line.startswith('#') ]

        header = lines.pop(0)
        columns = header.split()
        matrix = dict()

        for row in lines:
            entries = row.split()
            row_name = entries.pop(0)
            matrix[row_name] = dict()

            if len(entries) != len(columns):
                raise Exception('columns and rows counts does not match\n',
                                'file: {}'.format(self.matrix_path))

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
                raise KeyError('Bad key pair in substitution matrix: {}'
                               .format(key))
        else:
            return self.matrix[key]

    def get_name(self):
        """Returns substitution matrix path"""
        return self.name

    @staticmethod
    def matrices(matrix_name):
        """Default stored matrices"""
        if matrix_name == 'blossum62':
            return '''
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
            '''

        if matrix_name == 'nucleicmatrix':
            return '''
   A  C  G  T  U  X
A  2 -1 -1 -1 -1 -1
C -1  2 -1 -1 -1 -1
G -1 -1  2 -1 -1 -1
T -1 -1 -1  2 -1 -1
U -1 -1 -1 -1  2 -1
X -1 -1 -1 -1 -1  1
            '''

        return None



def select(select_list, target, sele, method):
    """Creates pymol selection object"""


    select_name = string.Formatter().vformat(sele, (), 
        SafeDict(
            method=method,
            target=re.sub(r'[^\w]', '', target)[:10],
            id=new_id(sele)))
    
    # empty select
    cmd.select(select_name, None)

    for select_tuple in select_list:
        model = select_tuple[0]
        chain = select_tuple[1]
        resi = select_tuple[2]

        # select /model/?/chain/resi-resi
        select_query = " | /{0}//{1}/{2}".format(model, chain, resi)

        # Execute and append selection to select_id
        cmd.select(select_name, select_name + select_query)


def new_id(sele):
    """Returns an incremented id if id token is requested"""
    if re.search(r'{id}', sele):
        stored.id += 1
    
    return stored.id


class SafeDict(dict):
    def __missing__(self, key):
        return key


def parse_str_to_list(string):
    """ Parse a string and return a list of values """

    # split by any separator , . / etc. excluding all letters, numbers and _
    raw_str = re.sub('([A-Za-z0-9_]+)', r'\1', string)

    # remove symbols '[', ']', ''', '"', white space
    raw_str = re.sub('[\[\]\'\"\s]', '', raw_str)

    return raw_str.split(',')

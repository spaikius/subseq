from pymol import cmd, stored
import re
import os


class BadParameterError(Exception):
    """ Bad Parameter Exception """
    pass


class BadRegExSyntaxError(Exception):
    """ Regular Expresion syntax Exception """
    pass


class InvalidMatrixFormatError(Exception):
    """ Invalid substitution matrix Exception """
    pass


def subseq_re(target='', chains='all', search_for='aminoAcids', first_only='False', models='all'):
    try:
        models, chains, search_for, first_only = check_parameters(target, chains, search_for, first_only, models)
    except BadParameterError:
        print("[Info]  Errors were found. Please see above messages for more information")
        return

    # Print manual if target value is [-h, --h, -he, --he, -hel, --hel, -help, --help]
    if re.match(r'^-{1,2}h(?:elp|el|e|)$', target, re.I):
        print('help')
        return

    data = Data(models, chains, search_for)
    search_results = None
    try:
        search_results = subseq_re_search(target.upper(), data, first_only, search_for)
    except BadRegExSyntaxError as msg:
        print("[Error] {0}".format(str(msg)))

    if search_results is not None:
        select(search_results, method='re', target=target)
    else:
        print("[Info] Nothing can be found for given target: {0}".format(target))


def subseq_la(target='', sub_matrix='', chains='all', search_for='aminoacids', first_only='False', gap_cost='10',
              min_score='51', models='all'):
    try:
        models, chains, search_for, first_only = check_parameters(target, chains, search_for, first_only, models, 
                                                                  sub_matrix, gap_cost, min_score)
    except BadParameterError:
        print("[Info]  Errors were found. Please see above messages for more information")
        return

    if re.match(r'^-{1,2}h(?:elp|el|e|)$', target, re.I):
        print('help')
        return

    data = Data(models, chains, search_for)

    try:
        search_results = subseq_la_search(target, data, sub_matrix, gap_cost, min_score, first_only)
    except KeyError as msg:
        print("[Error] {0}".format(str(msg)))
        return

    if search_results is not None:
        select(search_results, method='la', target=target)
    else:
        print("[Info] Nothing can be found for given target: {0}".format(target))


cmd.extend('subseq.re', subseq_re)
cmd.extend('subseq.la', subseq_la)
stored.selection_id = 1


# -----------------------------------------------------------------------------
# CLASS DATA

class Data:
    """
    This class is designed to extract data from pymol cmd.iterate command
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

    def __init__(self, models, chains, search_for):
        self.data = None
        self.models = models
        self.chains = chains
        self.search_for = search_for.lower()

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
        Initializes data structure self.data
        structure schema:
            model: {chain: {sequence: '', ids: list()}}

        model and chain names are in self.models (list) and
        self.chains (list) respectively
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
                main_atoms[0] = 'X'

    def replace_to_na_one_letter(self, atoms_dict):
        """Replaces all DNA or RNA nucleic acid code to 1 letter aa code"""
        for main_atoms in atoms_dict['main_atoms']:
            try:
                main_atoms[0] = self.na_one_letter[main_atoms[0]]
            except KeyError:
                main_atoms[0] = 'X'

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
        """Remove all blank attributes in self.data dictionary"""
        for model in self.data.keys():
            for chain in self.data[model].keys():
                # remove chain if sequence is empty
                if not self.data[model][chain]['sequence']:
                    del self.data[model][chain]
                # remove model if does not contain any chain
                if not self.data[model]:
                    del self.data[model]


# END OF CLASS DATA
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# REGULAR EXPRESSION SEARCH

def subseq_re_search(target, data, firstOnly, search_for):
    """
    work flow:
        1) create a RegExp object from target (function parameter)
        2) scan data (function parameter) by using RegExp object
        3) append information about match to match_list
        4) return match_list if its length is not 0 else return None
    """
    match_list = list()

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
    except:
        raise BadRegExSyntaxError('Bad syntax - target(RegExp): ' + str(target))

    # scan data by using RegExp object
    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain]['sequence']):

                start_pos = match.start()
                for _ in range(0, len(match.group())):
                    resi = data[model][chain]['ids'][start_pos]

                    match_list.append((model, chain, resi))

                    start_pos += 1

                if firstOnly:
                    break
            else:
                continue
            break

    # return match list if its length is not 0 else return None
    return match_list if len(match_list) != 0 else None


# END OF REGULAR EXPRESSION SEARCH
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LOCAL ALIGNMENT

def subseq_la_search(target, data, matrix, gap_cost, min_score, first_only):
    sub_matrix = SubMatrix(matrix)
    max_score = calculate_max_score(target, sub_matrix)
    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            sw = SmithWaterman(target, sequence, gap_cost, sub_matrix)

            # Skip if alignment best score is less than minimum passing score
            # print SW.get_best_score(), max_score, min_score
            if (float(sw.get_best_score()) / max_score) * 100 < float(min_score):
                continue

            for i, j in sw.get_coordinates():
                aligned_target, aligned_sequence, start_i, start_j = sw.get_traceback(i, j)

                start_pos = start_j - 1
                for _ in range(0, len(aligned_sequence.replace('-', ''))):
                    resi = data[model][chain]['ids'][start_pos]
                    match_list.append((model, chain, resi))
                    start_pos += 1

                alignment_string, identities, gaps, mismatches = \
                    create_alignment_string(aligned_target, aligned_sequence)

                print_alignment(model, chain, target, sequence, sub_matrix.get_name(), gap_cost
                                , sw.get_best_score(), max_score, identities, mismatches, gaps
                                , aligned_target, aligned_sequence, alignment_string
                                , start_i, start_j, data[model][chain]['ids'])
                
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


def print_alignment(model, chain, target, sequence, substitution_matrix_name, gap_cost
                    , alignment_score, max_score, identities, mismatches, gaps
                    , aligned_target, aligned_sequence, alignment_string
                    , target_start, subject_start, ids_list):
    """
    Prints BLAST like alignment for both sequences
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
# CLASS SmithWaterman

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
        self.score_matrix = [[0 for _ in range(len(sequence) + 1)] for _ in range(len(target) + 1)]

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
        """
        Looks for the next move during traceback.
        Moves are determined by the score of three upper-left, left and up
        in the self.score_matrix elements
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]
        achieved_score = self.score_matrix[i][j]
        diagonal = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if achieved_score == diagonal + int(self.sub_matrix[aa1, aa2]):
            # return diagonal move if diagonal move is greater than 0 else return END
            return 1 if diagonal > 0 else 0

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

        diagonal_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(0, diagonal_score, up_score, left_score)


# END OF CLASS SmithWaterman
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CLASS SubMatrix

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
                raise KeyError("Bad key pair in substitution matrix: {}".format(key))
        else:
            return self.matrix[key]

    def get_name(self):
        return self.name


# END OF CLASS SubMatrix
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# SELECT

def select(select_list, method, target):
    # select ID
    select_id = 'ss_' + str(method) + '_' + str(stored.selection_id) + '_' + str(target[:8])
    stored.selection_id += 1
    # empty select
    cmd.select(select_id, None)

    for select_tuple in select_list:
        model = select_tuple[0]
        chain = select_tuple[1]
        resi = select_tuple[2]

        # select /model/?/chain/resi-resi
        select_query = " | /{0}//{1}/{2}".format(model, chain, resi)

        # Execute and append selection to select_id
        cmd.select(select_id, select_id + select_query)


# END OF SELECT
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
def check_parameters(target, chains, search_for, first_only, models, sub_matrix=None, gap_cost=None, min_score=None):
    errors_found = False  # False - No errors were found, True - Errors were found

    all_models = cmd.get_names()
    # all_chains = get_all_chains(all_models)
    all_chains = list()

    for model in all_models:
        for chain in cmd.get_chains(model):
            all_chains.append(chain)

    # Remove all duplicates
    all_chains = list(set(all_chains))

    if target == '':
        print("[Error] parameter 'target' is not specified.")
        errors_found = True

    if models.lower() == 'all':
        models = all_models
    else:
        models = parse_str_to_list(models)

    if chains.lower() == 'all':
        chains = all_chains
    else:
        chains = [chain.upper() for chain in parse_str_to_list(chains)]

    for model in models:
        if model not in all_models:
            print("[Error] Model '{0}' does not exist.".format(model))
            errors_found = True

    for chain in chains:
        if chain not in all_chains:
            print("[Error] Chain '{0}' does not exist in any provided models".format(chain))
            errors_found = True

    if first_only in ['True', 'true', '1']:
        first_only = True
    elif first_only in ['False', 'false', '0']:
        first_only = False
    else:
        print("[Error] The 'firstOnly' parameter was not True, 1, False or 0.")
        errors_found = True

    if search_for.lower() in ['aminoacids', 'amino', 'aa']:
        search_for = 'aminoacids'
    elif search_for.lower() in ['nucleicacids', 'nucleic', 'na']:
        search_for = 'nucleicacids'
    else:
        print("[Error] The 'searchFor' parameter was not aminoAcid, amino, aa or nucleicAcid, nucleic, na")
        errors_found = True

    if sub_matrix:
        if sub_matrix == '':
            print("[Error] no substitution matrix provided")
            errors_found = True
        elif not os.path.isfile(sub_matrix):
            print("[Error] file '{0}' does not exist.".format(os.path.basename(sub_matrix)))
            errors_found = True

    if gap_cost:
        try:
            float(gap_cost)
        except TypeError:
            print("Gap cost must be type of int or float. Got {}".format(gap_cost))
            errors_found = True

    if min_score:
        try:
            if float(min_score) < 0 or float(min_score) > 100:
                print('minScore value must be between 0 and 100')

        except TypeError:
            print('minScore should be type of int or float')

    if errors_found:
        raise BadParameterError

    return models, chains, search_for, first_only


def get_all_chains(model):
    """ return a list of chains from specific model """

    chain_list = list()
    for chain in cmd.get_chains(model):
        chain_list.append(chain.upper())

    return chain_list


def parse_str_to_list(string):
    """ Parse a string and return a list of values """

    # split by any separator (, . / etc.) excluding all letters, numbers and _
    raw_str = re.sub('([A-Za-z0-9_]+)', r'\1', string)

    # remove symbols '[', ']', ''', '"', white space
    raw_str = re.sub('[\[\]\'\"\s]', '', raw_str)

    return raw_str.split(',')

# END OF HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# @module ss_LA_SEARCH.py
# @public functions: -
# @public variables: -
# @private functions: -
# @private variables: -


class SubMatrix:
    def __init__(self, matrix_path):
        self._load_matrix(matrix_path)

    def _load_matrix(self, matrix_path):
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
                raise Exception('Improper entry number in row')
            for column_name in columns:
                matrix[row_name][column_name] = entries.pop(0)

        self.matrix = matrix

    def get_score(self, a, b):
        if a not in self.matrix or b not in self.matrix[a]:
            raise Exception('Bad pair in substitution matrix: [%s, %s]' % (a, b))
        return self.matrix[a][b]


class SWA:
    def __init__(self, seq1, seq2, gap_cost, extend_cost):
        self.seq1 = [i for i in seq1]
        self.seq2 = [i for i in seq2]
        self.rows = len(seq1) + 1
        self.cols = len(seq2) + 1
        self.score_matrix = \
            [
                [0 for col in range(self.rows)]
                for row in range(self.cols)
            ]
        self.gap = gap_cost
        self.extend = extend_cost
        self.max_score = 0
        self.max_score_cord = None
        self.aligned_seq1 = None
        self.aligned_seq2 = None
        self.alignment = None

    def fill_data(self, sub_matrix):
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                score = self._calc_score(i, j, sub_matrix)
                self.score_matrix[i][j] = score

                if score > self.max_score:
                    self.max_score = score
                    self.max_score_cord = (i, j)

    def _calc_score(self, i, j, sub_matrix):
        aa1 = self.seq1[i - 1]
        aa2 = self.seq2[j - 1]

        similarity = sub_matrix.get_score(aa1, aa2)

        diag_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap
        left_score = self.score_matrix[i][j - 1] - self.gap

        return max(diag_score, up_score, left_score, 0)

    def traceback(self):
        a_seq1_list = list()
        a_seq2_list = list()
        i, j = self.max_score_cord
        move = self._next_move(i, j)

        while move != 0:
            # diag
            if move == 1:
                a_seq1_list.append(self.seq1[i - 1])
                a_seq2_list.append(self.seq2[j - 1])
                i -= 1
                j -= 1
            # up 
            if move == 2:
                a_seq1_list.append(self.seq1[i - 1])
                a_seq2_list.append('-')
                i -= 1
            # left
            if move == 3:
                a_seq1_list.append('-')
                a_seq2_list.append(self.seq2[j - 1])
                j -= 1

            move = self._next_move(i, j)

        a_seq1_list.append(self.seq1[i - 1])
        a_seq2_list.append(self.seq2[j - 1])

        self.aligned_seq1 = ''.join(reversed(a_seq1_list))
        self.aligned_seq2 = ''.join(reversed(a_seq2_list))

    def _next_move(self, i, j):
        diag = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if diag >= up and diag >= left:
            return 1 if diag != 0 else 0
        if up > diag and up > left:
            return 2 if up != 0 else 0
        if left > diag and left >= up:
            return 3 if left != 0 else 0

    def _alignment_seq(self):
        pass


def subseq_la(_target, _data, _matrix, _gap_cost, _extend_cost):
    sub_matrix = SubMatrix(_matrix)

    for model in _data.keys():
        for chain in _data[model].keys():
            alignment = SWA(_target, chain, _gap_cost, _extend_cost)
            alignment.fill_data(sub_matrix)
            alignment.traceback()

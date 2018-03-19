# @module ss_LA_SEARCH.py
# @public functions: -
# @public variables: -
# @private functions: -
# @private variables: -

import ss_SubMatrix as SubMatrix

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
        _alignment_seq()


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
        for aa1, aa2 in zip(self.aligned_seq1, self.aligned_seq2):
            if aa1 == aa2:
                self.alignment += '|'
            elif '-' in (aa1, aa2):
                self.alignment += ' '
            else:
                self.alignment += ':'


    def print_r(self):
        for i in range(0, alength, 60):
            seq1_slice = self.aligned_seq1[i:i+60]
            print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            print('             {0}'.format(self.alignment[i:i+60]))
            seq2_slice = self.aligned_seq2[i:i+60]
            print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            print()


def subseq_la(_target, _data, _matrix, _gap_cost, _extend_cost):
    sub_matrix = SubMatrix.SubMatrix(_matrix)

    for model in _data.keys():
        for chain in _data[model].keys():
            alignment = SWA(_target, chain, _gap_cost, _extend_cost)
            alignment.fill_data(sub_matrix)
            alignment.traceback()
            alignment.print_r()

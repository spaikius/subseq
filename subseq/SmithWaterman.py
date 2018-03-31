
class SmithWaterman:
    def __init__(self, target, sequence, gap_cost, matrix):
        self.target = [i for i in target]
        self.sequence = [i for i in sequence]
        self.gap_cost = gap_cost
        self.sub_matrix = matrix
        self.score_matrix = None

        self.create_score_matrix()
        self.max_score, self.max_score_cord = self.fill_matrix()
        self.find_traceback()
        self.construct_aligment_string()

    def create_score_matrix(self):
        rows = len(self.target) + 1
        cols = len(self.sequence) + 1

        self.score_matrix = [[0 for col in range(cols)] for row in range(rows)]

    def fill_matrix(self):
        max_score = None
        max_score_cord = (None, None)
        for i in range(1, len(self.score_matrix)):
            for j in range(1, len(self.score_matrix[i])):
                score = self.calc_score(i, j)

                self.score_matrix[i][j] = score

                if(score > max_score):
                    max_score = score
                    max_score_cord = (i, j)

        return max_score, max_score_cord

    def calc_score(self, i, j):
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]

        similarity = self.sub_matrix[aa1, aa2]

        diag_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(0, diag_score, up_score, left_score) 

    def find_traceback(self):
        aligned_seq1 = list()
        aligned_seq2 = list()

        end, diag, up, left = range(4)
        i, j = self.max_score_cord

        move = self.next_move(i, j)

        while move != end:
            if move == diag:
                aligned_seq1.append(self.target[i - 1])
                aligned_seq2.append(self.sequence[j - 1])

                i -= 1
                j -= 1

            elif move == up:
                aligned_seq1.append(self.target[i - 1])
                aligned_seq2.append('-')

                i -= 1

            elif move == left:
                aligned_seq1.append('-')
                aligned_seq2.append(self.sequence[j - 1])

                j -=1

            move = self.next_move(i, j)

        aligned_seq1.append(self.target[i - 1])
        aligned_seq2.append(self.sequence[j - 1])

        self.aligned_seq1 = ''.join(reversed(aligned_seq1))
        self.aligned_seq2 = ''.join(reversed(aligned_seq2))



    def next_move(self, i, j):
        diag = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if diag >= up and diag >= left:
            return 1 if diag != 0 else 0

        if up > diag and up > left:
            return 2 if up != 0 else 0

        if left > diag and left >= up:
            return 3 if left != 0 else 0

        raise Exception('Failed to find next move during find_traceback in {}'
                        .format(__name__))

    def construct_aligment_string(self):
        indentities, gaps, mismatches = 0, 0, 0
        aligment_string = ''

        for aa1, aa2 in zip(self.aligned_seq1, self.aligned_seq2):
            if aa1 == aa2:
                aligment_string += '|'
                indentities += 1

            elif '-' in (aa1, aa2):
                aligment_string += ' '
                gaps += 1

            else:
                aligment_string += ':'
                mismatches += 1

        self.indentities = indentities
        self.gaps = gaps
        self.mismatches = mismatches
        self.aligment_string = aligment_string

    def print_aligment(self):
        a_len = len(self.aligment_string)
        print("\nIndentities: {0}/{1} ({2:.1%})"
              .format(self.indentities, a_len, self.indentities / a_len))

        print("Mismatches: {0}/{1} ({2:.1%})"
              .format(self.mismatches, a_len, self.mismatches / a_len))

        print("Gaps: {0}/{1} ({2:.1%})"
              .format(self.gaps, a_len, self.gaps / a_len))

        for i in range(0, a_len, 60):
            seq1_slice = self.aligned_seq1[i : i + 60]
            seq2_slice = self.aligned_seq2[i : i + 60]
            align_slice = self.aligment_string[i : i + 60]

            print("Target  {0:<4} {1} {2:<4}".format(i+1, seq1_slice, i + len(seq1_slice)))
            print("             {0}".format(align_slice))
            print("Subject {0:<4} {1} {2:<4}".format(i+1, seq2_slice, i + len(seq2_slice)))

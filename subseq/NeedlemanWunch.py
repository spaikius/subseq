

class NeedlemanWunsch:
    """
    This class performs nucleotide or protein sequence alignment using
    the Needleman-Wunsch algorithm
    """
    def __init__(self, target, sequence, gap_cost, sub_matrix):
        self.target = target
        self.sequence = sequence
        self.gap_cost = float(gap_cost)
        self.sub_matrix = sub_matrix
        self.score_matrix = [[0 for _ in range(len(sequence) + 1)]
                             for _ in range(len(target) + 1)]
        self.init_score_matrix()
        self.fill_score_matrix()

    def init_score_matrix(self):
        """Initialization of scoring matrix"""
        for i in range(len(self.target) + 1):
            self.score_matrix[i][0] = -self.gap_cost * i
        for j in range(len(self.sequence) + 1):
            self.score_matrix[0][j] = -self.gap_cost * j

    def fill_score_matrix(self):
        """Fills core matrix with scores representing trial alignments
        of the two sequences
        """
        for i in range(1, len(self.score_matrix)):
            for j in range(1, len(self.score_matrix[i])):
                score = self.calculate_score(i, j)

                self.score_matrix[i][j] = score

    def calculate_score(self, i, j):
        """Calculates score for given i and j position in the score matrix
        The score is based on the upper-left, left and up elements
        """
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]

        similarity = float(self.sub_matrix[aa1, aa2])

        diagonal_score = self.score_matrix[i - 1][j - 1] + similarity
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(diagonal_score, up_score, left_score)

    def get_traceback(self):
        """Finds the optimal path through the score matrix.
        Returns constructed alignment strings for target and subject and
        values of i, j where alignment begins
        """
        aligned_target = list()
        aligned_subject = list()

        end, diagonal, up, left = range(4)

        i, j = len(self.target), len(self.sequence)

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

        if i == 0 or j == 0:
            # return END
            return 0

        if achieved_score == diagonal + int(self.sub_matrix[aa1, aa2]):
            # return diagonal move
            return 1

        if achieved_score == up - self.gap_cost:
            # return up move
            return 2

        if achieved_score == left - self.gap_cost:
            # return left move
            return 3

    def get_alignment_score(self):
        """Returns aligment score"""
        return self.score_matrix[-1][-1]

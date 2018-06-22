"""Description

This module provides a class for optimal local alignment using Smith-Waterman algorithm
Wikipedia link: https://en.wikipedia.org/wiki/Smith-Waterman_algortihm

"""
class SmithWaterman:
    """
    This class performs nucleotide or protein sequence (depending on given
    substitution matrix) alignment using the Smith-waterman algorithm
    """

    def __init__(self, target, sequence, gap_cost, sub_matrix):
        self.target = target
        self.sequence = sequence
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
        self.score_matrix = [[0 for _ in range(len(sequence) + 1)]
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

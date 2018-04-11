"""Description

This module provides a class for optimal local alignment using Smith-Waterman algorithm
Wikipedia link: https://en.wikipedia.org/wiki/Smith-Waterman_algortihm

"""
from Exceptions import InvalidPairError


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

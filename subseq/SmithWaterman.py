"""Description

This module provides a class for optimal local alignment using Smith-Waterman algorithm
Wikipedia link: https://en.wikipedia.org/wiki/Smith-Waterman_algortihm

SmithWaterman (class):
    Class description:
        This class performs nucleotide or protein sequence (depending on given 
        substitution matrix) alignment using the Smith-waterman algorithm

    Class attributes:
        __init__(target: str, sequence: str, gap_cost: float, matrix: SubMatrix, minscore: float):
            Constructor. Calls create_score_matrix() and fill_matrix()

        get_alignment_coordinates() -> list:
            return a list of starting and ending positions of aligned sequences

        calc_max_score() -> int:
            calculate the best possible alignemnt score for given target

        create_score_matrix():
            create a 2d matrix amd set all scores to 0
            
        fill_matrix():
            Fill self.score_matrix with scores representing trial alignments of the two sequences

        calc_score(i: int, j: int) -> int:
            Calculate score for given i and j position in the self.score_matrix
            The score is based on the upper-left, left and up elements in self.score_matrix

        find_traceback(i: int, j: int) -> str, str, int, int:
            For the given best score coordinate find the optimal path through the self.score_matrix.
            Return constructed alignment strings for both sequences and coordinate where alignment begins
            
        next_move(i: int, j: int) -> int
            Looks for the next move during traceback.
            Moves are determined by the score of three upper-left, left and up 
            in the self.score_matrix elements

        construct_alignment_string(aligned_seq1: str, aligned_seq2: str) -> str, int, int, int:
            Construct alignment string for both aligned sequences

            KTGTA
            :| :|  <-- alignment string
            PT-KA

            where ':' - mismatch, ' ' - gap, '|' - match
            return alignment string, match count, mismatch count, gap count
    
        print_data():
            Calls find_traceback(), construct_aligment_string(), print_aligment() for each best
            acumamulative score in self.score_matrix

        print_alignment():
            Prints BLAST like alginment for both sequences
"""


class SmithWaterman:
    def __init__(self, target, sequence, gap_cost, matrix, minscore):
        self.target = [i for i in target]
        self.sequence = [i for i in sequence]
        self.gap_cost = gap_cost
        self.sub_matrix = matrix
        self.minscore = minscore

        self.alignment_coordinates = list()
        self.score_matrix = None

        self.max_score = self.calc_max_score()

        self.create_score_matrix()
        self.best_score, self.best_score_coord_list = self.fill_matrix()

    def get_alignment_coordinates(self):
        for i, j in self.best_score_coord_list:
            end, diag, up, left = range(4)

            move = self.next_move(i, j)

            start_coord = None
            end_coord = j - 1

            while move != end:
                if move == diag:
                    i -= 1
                    j -= 1
                elif move == up:
                    i -= 1
                elif move == left:
                    j -= 1

                move = self.next_move(i, j)

            start_coord = j - 1

            self.alignment_coordinates.append((start_coord, end_coord))

        return self.alignment_coordinates

    def calc_max_score(self):
        max_score = 0
        for aa in self.target:
            max_score += int(self.sub_matrix[aa, aa])

        return max_score

    def create_score_matrix(self):
        rows = len(self.target) + 1
        cols = len(self.sequence) + 1

        self.score_matrix = [[0 for col in range(cols)] for row in range(rows)]

    def fill_matrix(self):
        best_score = None
        best_score_coord_list = list()
        for i in range(1, len(self.score_matrix)):
            for j in range(1, len(self.score_matrix[i])):
                score = self.calc_score(i, j)

                self.score_matrix[i][j] = score

                score_in_perc = (float(score) / self.max_score) * 100

                # skip if acummalitve score is lower than minimum score
                if self.minscore > score_in_perc:
                    continue

                if score == best_score:
                    best_score_coord_list.append((i, j))
                elif score > best_score:
                    best_score = score
                    best_score_coord_list = [(i, j)]

        return best_score, best_score_coord_list

    def calc_score(self, i, j):
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]

        similarity = self.sub_matrix[aa1, aa2]

        diag_score = self.score_matrix[i - 1][j - 1] + int(similarity)
        up_score = self.score_matrix[i - 1][j] - self.gap_cost
        left_score = self.score_matrix[i][j - 1] - self.gap_cost

        return max(0, diag_score, up_score, left_score)

    def find_traceback(self, i, j):
        aligned_seq1 = list()
        aligned_seq2 = list()

        end, diag, up, left = range(4)

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

                j -= 1

            move = self.next_move(i, j)

        aligned_seq1.append(self.target[i - 1])
        aligned_seq2.append(self.sequence[j - 1])

        aligned_seq1 = ''.join(reversed(aligned_seq1))
        aligned_seq2 = ''.join(reversed(aligned_seq2))

        return aligned_seq1, aligned_seq2, i, j

    def next_move(self, i, j):
        aa1 = self.target[i - 1]
        aa2 = self.sequence[j - 1]
        achieved_score = self.score_matrix[i][j]
        diag = self.score_matrix[i - 1][j - 1]
        up = self.score_matrix[i - 1][j]
        left = self.score_matrix[i][j - 1]

        if achieved_score == diag + int(self.sub_matrix[aa1, aa2]):
            return 1 if diag > 0 else 0

        if achieved_score == up - self.gap_cost:
            return 2 if up > 0 else 0

        if achieved_score == left - self.gap_cost:
            return 3 if left > 0 else 0

        raise Exception('Failed to find next move in {}'
                        .format(__name__))

    @staticmethod
    def construct_alignment_string(aligned_seq1, aligned_seq2):
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

    def print_data(self):
        for i, j in self.best_score_coord_list:

            aligned_seq1, aligned_seq2, seq1_start, seq2_start = \
                self.find_traceback(i, j)

            alignment_string, identities, gaps, mismatches = \
                self.construct_alignment_string(aligned_seq1, aligned_seq2)

            self.print_alignment(aligned_seq1
                                 , aligned_seq2
                                 , alignment_string
                                 , identities
                                 , gaps
                                 , mismatches
                                 , seq1_start
                                 , seq2_start)

    def print_alignment(self
                        , aligned_seq1
                        , aligned_seq2
                        , alignment_string
                        , identities
                        , gaps
                        , mismatches
                        , seq1_start
                        , seq2_start):

        a_len = len(alignment_string)

        print("\nScore:      {0}/{1} ({2:.1%})"
              .format(self.best_score, self.max_score, float(self.best_score) / self.max_score))

        print("Identities: {0}/{1} ({2:.1%})"
              .format(identities, a_len, float(identities) / a_len))

        print("Mismatches: {0}/{1} ({2:.1%})"
              .format(mismatches, a_len, float(mismatches) / a_len))

        print("Gaps:       {0}/{1} ({2:.1%})"
              .format(gaps, a_len, float(gaps) / a_len))

        for i in range(0, a_len, 60):
            seq1_slice = aligned_seq1[i: i + 60]
            seq2_slice = aligned_seq2[i: i + 60]
            align_slice = alignment_string[i: i + 60]
            target_start = i + seq1_start
            subject_start = i + seq2_start
            target_end = len([i for i in seq1_slice if i != '-']) + target_start - 1
            subject_end = len([i for i in seq1_slice if i != '-']) + subject_start - 1

            print("Target  {0:<4} {1} {2:<4}".format(target_start, seq1_slice, target_end))
            print("             {0}".format(align_slice))
            print("Subject {0:<4} {1} {2:<4}".format(subject_start, seq2_slice, subject_end))

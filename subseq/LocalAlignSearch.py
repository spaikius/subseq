"""Description
This module is designed to perform local alignment for each model chain sequence and given target
"""

import SubMatrix
import SmithWaterman


def subseq_la(target, data, matrix, gap_cost, min_score):
    """
    
    """
    sub_matrix = SubMatrix.SubMatrix(matrix)
    max_score = calculate_max_score(target, sub_matrix)
    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            SW = SmithWaterman.SmithWaterman(target, sequence, gap_cost, sub_matrix)

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

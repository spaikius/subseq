from __future__ import print_function

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
        target_start_index, subject_start_index, ids_list):
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
          .format(gaps, a_len, float(gaps) / a_len))

    for i in range(0, a_len, 60):
        target_slice = aligned_target[i: i + 60]
        subject_slice = aligned_sequence[i: i + 60]
        alignment_slice = alignment_string[i: i + 60]

        target_start = i + target_start_index
        subject_start = i + int(ids_list[subject_start_index - 1])

        target_end = len([i for i in target_slice if i != '-']) \
            + target_start - 1

        subject_end = len([i for i in subject_slice if i != '-']) \
            + subject_start - 1

        print("Target  {0:<4} {1} {2}"
              .format(target_start, target_slice, target_end))

        print(' ' * 13 + "{0}".format(alignment_slice))

        print("{0}/{1:<2} {2:<4} {3} {4}"
              .format(model, chain, subject_start, subject_slice, subject_end))

        print("\n")

    print('-' * 60)

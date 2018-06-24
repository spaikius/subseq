import logging

import alignment
import CallCounter
import subseq_parse
import subseq_select
import SubMatrix
import SmithWaterman
import Data


def subseq_local_alignment(
        targets, submatrix='blossum62', chains='all', search='aminoacids',
        firstonly='False', gapcost='10.', minscore='51.', models='all',
        sele='ss-{method}-{id}-{target}'):
    """
DESCRIPTION
    subseq.local - tool for searching target sequences using local alignment

USAGE
    subseq.local targets, [submatrix, [chains, [search, [firstonly, [gapcost,
                 minscore, [models, [sele]]]]]]]

IMPORTANT
    All modified amino or nucleic acids are replaced with: X

PARAMETERS
    targets=<list|FILE>     ; Target sequence
                              Examples:
                                - targets=KTGTAVU
                                - targets=SIS KATK AK
                                - targets=PATH/TO/TARGETS_FILE

    submatrix=<FILE>        ; Path to substitution matrix file
                              Default: blossum62

    chains=<list>           ; The list of chains
                              Examples:
                                - chains=A
                                - chains=A AT X Q
                              Default: all

    search=<str>            ; Search for nucleic acids or amino acids sequence
                                - for amino acids: aminoAcids, amino, aa
                                - for nucleic acids: nucleicAcids, nucleic, na
                              Default value: aminoAcids

    firstonly=<bool>        ; If firstonly is False (0) then select all matches
                              If firstonly is True  (1) then select first match
                              Default: False

    gapcost=<float>         ; The linear gap cost for local alignment
                              Default value: 10.00

    minscore=<float>        ; The minimum score in precentages for throwing off
                              low score alignments in alignment search
                              Example: minscore=75.25
                              Default value: 51.00 ( 51% )

    models=<list>           ; The list of models
                              Examples:
                                - models=5ara
                                - models=[5ara, 2cif, a4s2]
                              Default: all

    sele=<str>              ; Selection name
                              Tokens:
                                - {method} - used method for sequence search
                                - {target} - first 10 alpha-numeric symbols
                                - {id}     - id
                              Default: 'ss-{method}-{id}-{target}'

EXAMPLE
    subseq.local KTGT, blossum62, firstonly=True, gapcost=12.5, chains=A B

SEE ALSO
    subseq, subseq.global

SUBSEQ                          2018-06-01
    """
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.error = CallCounter.CallCounter(logging.error)

    targets = subseq_parse.parse_targets(targets)
    chains = subseq_parse.parse_chains(chains)
    search = subseq_parse.parse_search(search)
    firstonly = subseq_parse.parse_firstonly(firstonly)
    gapcost = subseq_parse.parse_gapcost(gapcost)
    minscore = subseq_parse.parse_minscore(minscore)
    models = subseq_parse.parse_models(models)

    if logging.error.counter is not 0:
        logging.info("{0} errors were found. ".format(logging.error.counter) +
                     "Please see above messages for more information")

        return

    if search is 'nucleicacids' and submatrix is 'blossum62':
        submatrix = 'nucleicmatrix'

    data = Data.Data(models, chains, search, replace_with='X')

    for target in targets:
        try:
            search_results = subseq_la_search(target, data, submatrix,
                                              gapcost, minscore, firstonly)
        except Exception as e:
            logging.error("{0}".format(e))
            continue

        if search_results is not None:
            subseq_select.select(search_results, target, sele, method='local')

        else:
            logging.info("Nothing can be found for given target: {0}"
                         .format(target))


def subseq_la_search(target, data, matrix, gap_cost, min_score, first_only):
    # Substitution matrix
    sub_matrix = SubMatrix.SubMatrix(matrix)

    # The maximum score for given target
    max_score = alignment.calculate_max_score(target, sub_matrix)

    match_list = list()

    for model in data.keys():
        for chain in data[model].keys():
            sequence = data[model][chain]['sequence']
            sw = SmithWaterman.SmithWaterman(target, sequence, gap_cost, sub_matrix)

            # Skip if alignment best score is less than minimum passing score
            if float(sw.get_best_score()) / max_score * 100 < min_score:
                continue

            for i, j in sw.get_coordinates():

                aligned_target, aligned_sequence, start_i, start_j =\
                    sw.get_traceback(i, j)

                start_pos = start_j - 1

                for _ in range(0, len(aligned_sequence.replace('-', ''))):
                    resi = data[model][chain]['ids'][start_pos]
                    match_list.append((model, chain, resi))
                    start_pos += 1

                alignment_string, identities, gaps, mismatches = \
                    alignment.create_alignment_string(aligned_target, aligned_sequence)

                alignment.print_alignment(
                    model, chain, target, sequence, sub_matrix.get_name(),
                    gap_cost, sw.get_best_score(), max_score, identities,
                    mismatches, gaps, aligned_target, aligned_sequence,
                    alignment_string, start_i, start_j,
                    data[model][chain]['ids'])

                if first_only:
                    break
            else:
                continue
            break

    return match_list if len(match_list) != 0 else None
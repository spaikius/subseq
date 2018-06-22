import logging
import re

import subseq_parse
import subseq_select
import CallCounter
import Data

def subseq_re(
        targets='', chains='all', search='aminoAcids', firstonly='False',
        models='all', sele='ss-{method}-{id}-{target}'):
    """
DESCRIPTION
    subseq - tool for searching target sequences using Regular Expressions

USAGE
    subseq targets, [chains, [search, [firstonly, [models, [sele]]]]]

IMPORTANT
    All modified amino or nucleic acids are replaced with: X

    If quantifier {n,m} is used in Regular Expressions then target value should
    be within parentheses

                    targets= (GT{3,}A{,4}(L|G){2,10})
                             ^~~~~~~~~~~~~~~~~~~~~~~^
PARAMETERS
    targets=<list|FILE>     ; Target sequence
                              Examples:
                                - targets=KTGTAVU
                                - targets=(TATA.{3,5}ATG(.{3,4}){3,})
                                - targets=SIS KATK (AK{3,4})
                                - targets=PATH/TO/TARGETS_FILE

    chains=<list>           ; The list of chains
                              Examples:
                                - chains=A
                                - chains=A AT X Q
                              Default: all

    search=<str>            ; Search for nucleic acids or amino acids sequence
                                - for amino acids: aminoacids, amino, aa
                                - for nucleic acids: nucleicacids, nucleic, na
                              Default value: aminoacids

    firstonly=<bool>        ; If firstonly is False (0) then select all matches
                              If firstonly is True  (1) then select first match
                              Default: False

    models=<list>           ; The list of models
                              Examples:
                                - models=5ara
                                - models=5ara 2cif a4s2
                              Default: all

    sele=<str>              ; Selection name
                              Tokens:
                                - {method} - used method for sequence search
                                - {target} - first 10 targets alphabet letters
                                - {id}     - id
                              Default: 'ss-{method}-{id}-{target}'

EXAMPLE
    subseq KTGT (KT{2,4}), A B C, firstonly=True, search=nucleicacids

SEE ALSO
    subseq.local, subseq.global

SUBSEQ                          2018-06-01
    """

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.error = CallCounter.CallCounter(logging.error)

    targets = subseq_parse.parse_targets(targets)
    chains = subseq_parse.parse_chains(chains)
    search = subseq_parse.parse_search(search)
    firstonly = subseq_parse.parse_firstonly(firstonly)
    models = subseq_parse.parse_models(models)

    if logging.error.counter is not 0:
        logging.info("{0} errors were found. ".format(logging.error.counter) +
                     "Please see above messages for more information")

        return

    data = Data.Data(models, chains, search, replace_with='X')

    search_results = None

    for target in targets:
        try:
            search_results = subseq_re_search(target, data, firstonly, search)

        except Exception as e:
            logging.warning("RegExp for {0}: {1}".format(target, e))
            continue

        if search_results is not None:
            subseq_select.select(search_results, target, sele, method='re')

        else:
            logging.info("Nothing can be found for given target: {0}"
                         .format(target))


def subseq_re_search(target, data, first_only, search_for):
    """
    work flow:
        1) create a RegExp object
        2) scan data by using RegExp object
        3) append information about match to match_list
        4) return match_list if its length is not 0 else return None
    """
    match_list = list()

    target = target.strip("'()\"")

    # Replace wildcards
    if search_for == 'nucleicacids':
        target = target.replace('R', '[AG]')
        target = target.replace('Y', '[CT]')
        target = target.replace('S', '[GC]')
        target = target.replace('W', '[AT]')
        target = target.replace('K', '[GT]')
        target = target.replace('M', '[AC]')
        target = target.replace('B', '[CGT]')
        target = target.replace('D', '[AGT]')
        target = target.replace('H', '[ACT]')
        target = target.replace('N', '.')

    # RegExp validation
    try:
        # re.I - ignore case sensitive
        re_target = re.compile(target, re.I)
    except Exception:
        raise

    # scan data by using RegExp object
    for model in data.keys():
        for chain in data[model].keys():
            for match in re_target.finditer(data[model][chain]['sequence']):
                start_pos = match.start()

                for _ in range(0, len(match.group())):
                    resi = data[model][chain]['ids'][start_pos]

                    match_list.append((model, chain, resi))

                    start_pos += 1

                if first_only:
                    break
            else:
                continue
            break

    return match_list if len(match_list) else None
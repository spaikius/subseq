from pymol import cmd

import os
import sys


def __init__(self):
    path = os.path.dirname(__file__)
    sys.path.append(path)

    import ParseKwargs
    import Data
    import RegExSearch
    import LocalAlignSearch
    import Select
    import Exceptions


def subseq(*argv, **kwargs):
    # If user just typed 'subseq'
    # kwargs contains `_self` key, so it's length is 1, not 0
    if len(argv) == 0 and len(kwargs) == 1:
        print("For help type: subseq help")
        return

    # Print usage manual if user typed `subseq help` 
    if 'help' in argv:
       print(usage_message)
       return

    # Print warrning message if user forgot to seperate kwargs with comma.
    # Or targets value is not within parentheses if RegExp quantifier `{n,m}` is used
    if len(kwargs) != 1 and len(argv) > 0:
        print(bad_arguments_message)
        return

    try:
        parameters = ParseKwargs.ParseKwargs(kwargs)

        data = Data.Data(parameters['models'], parameters['chains'])

        search_result = None
        if parameters['algorithm'] == 're':
            search_result = RegExSearch.subseq_re(parameters['target'], data)
        elif parameters['algorithm'] == 'la':
            search_result = LocalAlignSearch.subseq_la(parameters['target']
                                                       , data
                                                       , parameters['submatrix']
                                                       , parameters['gapcost']
                                                       , parameters['minscore'])

        if search_result is not None:
            Select.select(search_result)
        else:
            print("Nothing can be found for given target: {0}".format(parameters['target']))

    except (  Exceptions.BadParameterError
            , Exceptions.InvalidPairError
            , Exceptions.BadRegExSyntaxError
            , Exceptions.NoModelsError
            , Exceptions.InvalidMatrixFormatError
            , Exceptions.InvalidPairError ) as e:
        print(' [Error] ' + str(e))

    except Exception:
        # Unknown exception occured.
        # Reraise exception and let Pymol to handle it. 
        raise

cmd.extend('subseq', subseq)

usage_message = """
Usage: subseq target=<str>, algorithm=<str>, submatrix=<str>, models=<list>
              , chains=<list>, gapcost=<float>, minscore=<float>

Exaple usage: subseq target=KTGT, algorithm=la, chains=[A, B, T], submatrix=PATH/TO/MATRIX

!!! Important !!!
Please note: each keyword parameter should be seperated with comma (,)
Please note: target value should be within parentheses if quantifier {n,m} is used in Regular Expressions.
             Example: target=(GT{3,})
                             ^      ^


Parameters:
    help                            ; Prints usage manual

Keyword parameters:
    target=<str>        Required    ; Target sequence
                                      Examples:
                                       If algorithm type is re:
                                        - target=KTGTAVU
                                        - target=(TATA.{3,5}ATG(.{3,4}){3,})
                                       If algorithm type is la:
                                        - target=KTGAT


    algorithm=<str>     Optional    ; Algorithm search type.
                                      - 're' for Regular Expression
                                      - 'la' for local alignment. Smith-Waterman 
                                      Default value: 're'

    models=<list>       Optional    ; The list of models that will be used for target search.
                                      If the list is not provided, then search for target will be 
                                      performed in all available models
                                      Example: models=[5ara, 2cif, a4s2]

    chains=<list>       Optional    ; The list of chains that will be used for target search.
                                      If the list is not priveded, then search for target will be
                                      performed in all available model chains
                                      Example: chains=[A, T, X, Q]

    submatrix=<PATH>    Optional    ; Path to substitution matrix for local alignment
                                      Default substitution matrix: Blossum62

    gapcost=<float>     Optional    ; The linear gap cost for local alignment
                                      Default value: 10

    minscore=<float>    Optional    ; The minimum score in precentages for throwing off low score alignments
                                      in local alignment search.
                                      Default value: 51.00 ( 51% )
                                      Example: minscore=75.25
"""

bad_arguments_message = """
Found nonseperated keyword arguments or target= value is not within parentheses.

Please check if all keyword arguments are seperated with comma (,) or target= value
is within parentheses if Regular Expression quantifier is used.

For help type: subseq help
"""
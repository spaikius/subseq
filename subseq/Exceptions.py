"""
This module is for custom exceptions
"""

class InvalidPairError(Exception):
    """ Invalid key pair Exception """
    pass

class NoModelsError(Exception):
    """ No models found Exception """
    pass

class BadParameterError(Exception):
    """ Bad Parameter Exception """
    pass

class BadRegExSyntaxError(Exception):
    """ Regular Expresion syntax Exception """
    pass

class InvalidMatrixFormatError(Exception):
    """ Invalid substitution matrix Exception """
    pass

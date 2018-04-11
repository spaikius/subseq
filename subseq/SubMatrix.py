"""Description
This module provides a class for generating substitution martrix from a file
"""

import os
from Exceptions import InvalidMatrixFormatError, InvalidPairError


class SubMatrix:
    def __init__(self, matrix_path):
        self.matrix = None
        self.name = os.path.basename(matrix_path)
        self.load_matrix(matrix_path)

    def load_matrix(self, matrix_path):
        with open(matrix_path, 'r') as fh:
            matrix = fh.read()

        lines = matrix.strip().split('\n')
        # remove comments
        lines = [line for line in lines if line[0] != '#']

        header = lines.pop(0)
        columns = header.split()
        matrix = dict()

        for row in lines:
            entries = row.split()
            row_name = entries.pop(0)
            matrix[row_name] = dict()

            if len(entries) != len(columns):
                raise InvalidMatrixFormatError('columns and rows counts does not match\n file: {}'
                                .format(self.matrix))
            for column_name in columns:
                matrix[row_name][column_name] = entries.pop(0)

        self.matrix = matrix

    def __getitem__(self, key): 
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

    def get_name(self):
        return self.name

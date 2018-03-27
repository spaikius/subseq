"""Description

This module is designed to get a score from substitution matrix

Module attribute:
    SubMatrix (class):
        class attributes:
            - __init__(matrix_path: str): typed
                constructor
                
            - load_matrix(matrix_path: str): typed
                loads substitution matrix from given path in matrix_path
    
            - __getitem__(key: tuple/str): typed
                returns value for given key or keys set

        usage:
            matrix = SubMatrix(path_to_matrix)
            matrix['A']['L'] 
            # or  
            matrix['A', 'L']

"""


class SubMatrix:
    def __init__(self, matrix_path):
        self.matrix = None
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
                raise Exception('Improper entry number in row')
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
                 raise Exception("SubMatrix: Bad key {}".format(key))
        else:
            return self.matrix[key]

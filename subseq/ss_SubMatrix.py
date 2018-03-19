
class SubMatrix:
    def __init__(self, matrix_path):
        self._load_matrix(matrix_path)

    def _load_matrix(self, matrix_path):
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

    def get_score(self, a, b):
        if a not in self.matrix or b not in self.matrix[a]:
            raise Exception('Bad pair in substitution matrix: [%s, %s]' % (a, b))
        return self.matrix[a][b]
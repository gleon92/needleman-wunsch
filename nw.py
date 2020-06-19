import pandas as pd


def maximo(a, b, c):
    Max = a
    if b > Max:
        Max = b
    if c > Max:
        Max = c
        if b > c:
            Max = b
    return Max


def convert_csv(matrix_file):
    df = pd.read_csv(matrix_file, sep=',', header=None)
    data = df.values
    bases = ['A', 'C', 'G', 'T']
    matriz = {}
    for i, base_i in enumerate(bases):
        for j, base_j in enumerate(bases):
            matriz[base_i, base_j] = data[i][j]
    return matriz


def needleman_wunsch2(seq_1, seq_2, gap_penalty, matriz_sustitucion):

    def _traceback(t, r, str1, str2, x, y, s1='', s2=''):
        if x > 0 or y > 0:
            c = t[x][y]
            u = c == (t[x][y-1] + gap_penalty)
            l = c == (t[x-1][y] + gap_penalty)
            ul = c == (t[x - 1][y - 1] + matrix[str1[x-1], str2[y-1]])
            if ul:
                _traceback(t, r, str1, str2,
                           x - 1, y - 1, str1[x - 1] + s1, str2[y - 1] + s2)
            if l:
                _traceback(t, r, str1, str2,
                           x - 1, y, str1[x - 1] + s1, '-' + s2)
            if u:
                _traceback(t, r, str1, str2,
                           x, y - 1, '-' + s1, str2[y - 1] + s2)
        else:
            r.append((s1, s2))

    matrix = convert_csv(matriz_sustitucion)
    matriz_alineamiento = [[0 for j in range(len(seq_2) + 1)] for i in range(len(seq_1) + 1)]

    for i in range(len(seq_1)):
        matriz_alineamiento[0][i + 1] = (i + 1) * gap_penalty
    for j in range(len(seq_2)):
        matriz_alineamiento[j + 1][0] = (j + 1) * gap_penalty
    for i, val_i in enumerate(seq_1):
        for j, val_j in enumerate(seq_2):
            matriz_alineamiento[i + 1][j + 1] = maximo(matriz_alineamiento[i][j] + matrix[val_i, val_j],
                                                       matriz_alineamiento[i][j+1] + gap_penalty,
                                                       matriz_alineamiento[i+1][j] + gap_penalty)

    results = []
    _traceback(matriz_alineamiento, results, seq_1, seq_2, i+1, j+1)

    for i in matriz_alineamiento:
        print(i)

    for n, r in enumerate(results):
        print('Resultado N°: ' + str(n+1))
        print(r[1])
        print(r[0])


needleman_wunsch2('AGC', 'AAG', -5, 'matriz_sustitucion.csv')
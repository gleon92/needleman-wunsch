import pandas as pd
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo


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


def needleman_wunsch(seq_1, seq_2, gap_penalty, matriz_sustitucion=None):
    def backtrack(t, r, str1, str2, x, y, s1='', s2=''):
        if x > 0 or y > 0:
            c = t[x][y]
            arriba = c == (t[x][y - 1] + gap_penalty)
            izquierda = c == (t[x - 1][y] + gap_penalty)
            diagonal = c == (t[x - 1][y - 1] + (matrix[str1[x - 1], str2[y - 1]] if (str1[x - 1], str2[y - 1]) in matrix else matrix[str2[y - 1], str1[x - 1]]))
            if diagonal:
                backtrack(t, r, str1, str2,
                          x - 1, y - 1, str1[x - 1] + s1, str2[y - 1] + s2)
            if izquierda:
                backtrack(t, r, str1, str2,
                          x - 1, y, str1[x - 1] + s1, '-' + s2)
            if arriba:
                backtrack(t, r, str1, str2,
                          x, y - 1, '-' + s1, str2[y - 1] + s2)
        else:
            r.append((s1, s2))

    if matriz_sustitucion:
        matrix = convert_csv(matriz_sustitucion)
    else:
        matrix = MatrixInfo.blosum62
    score_matrix = [[0 for j in range(len(seq_2) + 1)] for i in range(len(seq_1) + 1)]

    for i in range(len(seq_1)):
        score_matrix[0][i + 1] = (i + 1) * gap_penalty
    for j in range(len(seq_2)):
        score_matrix[j + 1][0] = (j + 1) * gap_penalty
    for i, val_i in enumerate(seq_1):
        for j, val_j in enumerate(seq_2):
            score_matrix[i + 1][j + 1] = maximo(
                score_matrix[i][j] + (matrix[val_i, val_j] if ((val_i, val_j) in matrix) else matrix[val_j, val_i]),
                score_matrix[i][j + 1] + gap_penalty,
                score_matrix[i + 1][j] + gap_penalty)

    results = []
    score = score_matrix[i+1][j +1]
    backtrack(score_matrix, results, seq_1, seq_2, i + 1, j + 1)

    for i in score_matrix:
        print(i)

    for n, r in enumerate(results):
        print('Resultado N°: ' + str(n + 1))
        print(r[1])
        print(r[0])
        print('Score: ' + str(score))


in_handle = open('flna.fasta')
record_iterator = SeqIO.parse(in_handle, "fasta")
rec_1 = next(record_iterator).upper()
rec_2 = next(record_iterator).upper()
needleman_wunsch('AGC', 'AAG', -5, 'matriz_sustitucion.csv')
# needleman_wunsch(rec_1, rec_2, -5)

from Bio.SubsMat import MatrixInfo as mi
import numpy as np
from itertools import combinations


class traceback:
    """Class used for reconstruction of the optimal alignment path"""
    def __init__(self, point, MSA1, MSA2, partialMSA1, partialMSA2, matrix):
        self.point = point
        self.MSA1 = MSA1
        self.MSA2 = MSA2
        self.partialMSA1 = partialMSA1
        self.partialMSA2 = partialMSA2
        self.matrix = matrix
        self.output = []

    def step(self):

        if(self.point == (0, 0)):
            self.output = [x for x in self.partialMSA1] + [x for x in self.partialMSA2]
            return
        if(self.matrix[self.point[0]][self.point[1]] == 4
           or self.matrix[self.point[0]][self.point[1]] == 5
           or self.matrix[self.point[0]][self.point[1]] == 6
           or self.matrix[self.point[0]][self.point[1]] == 7):
            #self.partialMSA1 = self.seqA[self.point[0] - 1] + self.partialMSA1
            self.partialMSA1 = [x[self.point[0] - 1]+y for (x,y) in zip(self.MSA1, self.partialMSA1)]
            self.partialMSA2 = [x[self.point[1] - 1] + y for (x, y) in zip(self.MSA2, self.partialMSA2)]
            self.point = (self.point[0]-1, self.point[1] - 1)
            self.step()
        if (self.matrix[self.point[0]][self.point[1]] == 2 or self.matrix[self.point[0]][self.point[1]] == 3):
            #from the left, gap in MSA1
            #self.partialMSA1 = "-" + self.partialMSA1
            self.partialMSA1 = ["-" + x for x in self.partialMSA1]
            self.partialMSA2 = [x[self.point[1] - 1] + y for (x, y) in zip(self.MSA2, self.partialMSA2)]
            self.point = (self.point[0], self.point[1]-1)
            self.step()
        if (self.matrix[self.point[0]][self.point[1]] == 1):
            #from above, gap in B
            #self.partialMSA1 = self.MSA1[self.point[0] - 1] + self.partialMSA1
            self.partialMSA1 = [x[self.point[0] - 1] + y for (x, y) in zip(self.MSA1, self.partialMSA1)]
            #self.partialMSA2 = "-" + self.partialMSA2
            self.partialMSA2 = ["-" + x for x in self.partialMSA2]
            self.point = (self.point[0] -1, self.point[1])
            self.step()




def readValue(query):
    """Reads Blosum62 value of a pair of amino acids"""
    if query in mi.blosum62.keys():
        return mi.blosum62[query]
    else:
        return mi.blosum62[(query[1],query[0])]

def scoreColumns(MSA1, msa1_index, MSA2, msa2_index):
    """Returns the score of matching whole alignment columns, analogous to matching amino acids"""
    gap_penalty = -7
    score = 0
    column_1 = [row[msa1_index] for row in MSA1]
    # for row in MSA2:
    #     print(row, msa2_index)
    column_2 = [row[msa2_index] for row in MSA2]
    for residue in column_1:
        for residue2 in column_2:
            if residue == "-" or residue2 == "-":
                score = score + gap_penalty
            else:
                score = score + readValue((residue, residue2))

    return score


def NeedlemanWunschMSA(MSA1, MSA2):
    """Performs dynamic programming alignment of alignments"""
    gap_penalty = -7
    paths = np.zeros((len(MSA1[0]) + 1, len(MSA2[0]) + 1), int)
    for i in range(1, len(MSA2[0])):
        paths[0][i] = 2
    for i in range(1, len(MSA1[0])):
        paths[i][0] = 1

    H = np.zeros((len(MSA1[0]) + 1, len(MSA2[0]) + 1), int)
    for i in range(1, len(MSA2[0])+1):
        H[0][i] = i * gap_penalty
    for i in range(1, len(MSA1[0])+1):
        H[i][0] = i * gap_penalty



    for i in range(1, len(MSA1[0]) + 1):
        for j in range(1, len(MSA2[0]) + 1):
            #optimalValue = max(H[i-1][j-1] + readValue((MSA1[i - 1], MSA2[j - 1])), H[i - 1][j] + gap_penalty, H[i][j - 1] + gap_penalty)
            optimalValue = max(H[i - 1][j - 1] + scoreColumns(MSA1, i - 1, MSA2, j - 1), H[i - 1][j] + gap_penalty,
                               H[i][j - 1] + gap_penalty)

            H[i][j] = optimalValue
            #diagonal +4, right arrow +2, down arrow +1
            if optimalValue == H[i-1][j-1] + scoreColumns(MSA1, i - 1, MSA2, j - 1):
                #diagonal
                paths[i][j] = paths[i][j] + 4
            if optimalValue == H[i-1][j] + gap_penalty:
                #from above
                paths[i][j] = paths[i][j] + 1
            if optimalValue == H[i][j-1] + gap_penalty:
                #from the left
                paths[i][j] = paths[i][j] + 2
    print("Optimal value: ", H[len(MSA1[0])][len(MSA2[0])])
    traversal = traceback((len(MSA1[0]), len(MSA2[0])), MSA1, MSA2, ["" for seq in MSA1], ["" for seq in MSA2], paths)
    traversal.step()

    return traversal.output
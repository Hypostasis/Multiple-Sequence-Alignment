from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as mi
from itertools import combinations
import numpy as np

def readBlosum(query, scoring_matrix = mi.blosum62):
    """Returns the substitution score of a pair of proteins"""
    if query in scoring_matrix.keys():
        return scoring_matrix[query]
    else:
        return scoring_matrix[(query[1], query[0])]

def compute_sum_of_pairs(alignment, scoring_matrix = mi.blosum62, gap_penalty = -1):
    """Returns the sum of pairs evaluation score of the given alignment"""
    score = 0
    #print(alignment, alignment.shape)
    for i, column in enumerate(alignment.T):
        for pair in combinations(column, 2):
            if "-" in pair:
                if pair[0] == pair[1]:
                    #score of (-,-) is 0
                    continue
                score = score + gap_penalty
                continue
            score = score + readBlosum(pair, scoring_matrix)


    return score

def evaluate_clustal(filename, function = compute_sum_of_pairs):
    """reads a clustal file and prints its sum-of-pairs score to the console"""
    alignment = []
    for seq_record in SeqIO.parse(filename, "clustal"):
        alignment.append(seq_record)
    alignment_array = np.array(alignment)
    score = compute_sum_of_pairs(alignment_array)
    print(score)


evaluate_clustal("test_files/testcase1-output-clustal", compute_sum_of_pairs)
evaluate_clustal("test_files/testcase1-output-tcoffee", compute_sum_of_pairs)
evaluate_clustal("test_files/testcase2-output-clustal", compute_sum_of_pairs)
evaluate_clustal("test_files/testcase2-output-tcoffee", compute_sum_of_pairs)
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as mi
from itertools import combinations
import numpy as np
from random import random

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
import random

def random_AA_seq(length):
    return ''.join(random.choice('ACDEFGHIKLMNPQRSTVWY') for i in range(length))

def random_AA():
    return random.choice('ACDEFGHIKLMNPQRSTVWY')


def mutate(sequence, number_of_sequences, p_insert=0.049, p_del=0.078, p_sub=0.051):
    output = []
    for i in range(number_of_sequences):
        seq = ""
        for amino_acid in sequence:
            if random.random() <= p_del:
                continue
            if random.random() <= p_del + p_insert:
                # insertion
                seq = seq + amino_acid
                seq = seq + random_AA()
                continue
            if random.random() <= p_del + p_insert + p_sub:
                seq = seq + random_AA()
                continue
            seq = seq + amino_acid
        output.append(seq)

    with open('testcase.fasta', 'w') as file:
        for i, read in enumerate(output):
            file.write(">Sequence" + str(i) + "\n")
            file.write(read+"\n")

    return output

def generate_testcase(length, number_of_sequences):
    random_sequence = random_AA_seq(length)
    mutate(random_sequence, number_of_sequences)
    return

generate_testcase(100000, 10)

# evaluate_clustal("test_files/testcase1-output-clustal", compute_sum_of_pairs)
# evaluate_clustal("test_files/testcase1-output-tcoffee", compute_sum_of_pairs)
# evaluate_clustal("test_files/testcase2-output-clustal", compute_sum_of_pairs)
# evaluate_clustal("test_files/testcase2-output-tcoffee", compute_sum_of_pairs)
# evaluate_clustal("test_files/testcase3-output-clustal", compute_sum_of_pairs)
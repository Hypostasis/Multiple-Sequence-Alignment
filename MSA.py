from tkinter import Tk, Label, Button, Entry, Checkbutton, IntVar
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo as mi
from Bio.Align import AlignInfo
import numpy as np
import tests
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import combinations
from NeedlemanWunschMSA import NeedlemanWunschMSA

def levenshtein(s1, s2):
    """computes Levenshtein distance"""
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[
                             j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def distance(seq1, seq2):
    """Computes distance between two sequences"""
    return levenshtein(seq1, seq2)


def save_msa_to_file(msa, filename = "output.txt"):

    msa_array = np.array(msa)
    print("Computing sum of pairs...")
    score = tests.compute_sum_of_pairs(msa_array, mi.blosum62)
    print(score)
    print("Saving to file...")
    file = open(filename, "w")
    for x in msa:
        file.write(str(x))
        file.write("\n")

    for i in range(len(msa_array[0])):
        conservative = True
        for j in range(len(msa_array)):
            if msa_array[j][i] != msa_array[0][i]:
                conservative = False
        if conservative:
            file.write("*")
        else:
            file.write(" ")
    file.write("\n")
    file.write("Sum-of-pairs score: " + str(score))
    file.write("\n")
    file.close()
    print("File saved.")
    return score

def merge(left, right):
    #left and right are lists of sequences representing a partial MSA
    print("TO MERGE:")
    print(left)
    print(right)
    merged = [a+b for a, b in zip(left, right)]
    print("MERGED:")
    print(merged)
    return merged


def dca(sequences, l_min, match_score, mismatch_penalty, gap_penalty, extension_penalty):
    divide = False
    for x in sequences:
        if len(x) > l_min:
            divide = True

    if divide:
        left = [x[:len(x) // 2] for x in sequences]
        right = [x[len(x) // 2:] for x in sequences]
        print("LEFT", left)
        print("RIGHT", right)
        left_dca = dca(left, l_min, match_score, mismatch_penalty, gap_penalty, extension_penalty)
        right_dca = dca(right, l_min, match_score, mismatch_penalty, gap_penalty, extension_penalty)
        return merge(left_dca, right_dca)
    else:
        print("TO ALIGN:")
        print(sequences[0], sequences[1])
        alignment = pairwise2.align.globalms(sequences[0], sequences[1], match_score, mismatch_penalty, gap_penalty,
                                             extension_penalty, one_alignment_only=True)[0]
        alignment = [alignment[0], alignment[1]]
        print("PARTIAL ALIGNMENT:")
        print(alignment)
        return alignment

class Node:

    def __init__(self, msa):
        self.msa = msa
        self.compute_consensus()

    def compute_consensus(self):
        align = MultipleSeqAlignment(Gapped(IUPAC.extended_protein, "-"))
        for i, seq in enumerate(self.msa):
            align.add_sequence(str(i), str(seq))
        summary_align = AlignInfo.SummaryInfo(align)
        self.consensus = summary_align.gap_consensus(threshold=0, ambiguous="-")

    def __repr__(self):
        repr = "MSA:" + str(self.msa) + "CONSENSUS:" + str(self.consensus)
        return repr


def merge_nodes(node1, node2):
    """this function should perform optimal alignment of alignments (e.g. modified Needleman-Wunsch)"""
    msa1 = node1.msa
    msa2 = node2.msa
    print("Merging nodes...")
    if len(msa2[0]) > len(msa1[0]):
        msa1, msa2 = msa2, msa1

    if len(msa1) == 1 and len(msa2) == 1:
        alignment = pairwise2.align.globalms(msa1[0], msa2[0], 1, -1, -1, -1, one_alignment_only=True)[0]
        print(alignment)
        new_msa = [str(alignment[0]), str(alignment[1])]
        return Node(new_msa)

    new_msa = NeedlemanWunschMSA(msa1, msa2)
    return Node(new_msa)

class MSA:
    #the variable sequences is a list of tuples (sequence_description, sequence)
    sequences = []
    def __init__(self, master):
        self.master = master
        master.title("Multiple Sequence Alignment")

        self.label = Label(master, text = "Enter your sequences here or as a FASTA file:")
        self.label.pack()
        self.field_input_file_name = Entry(text="input")
        self.field_input_file_name.pack()
        self.align_button = Button(master, text="Load FASTA", command=self.load_fasta)
        self.align_button.pack()
        self.field_S1 = Entry()
        self.field_S1.pack()
        self.field_S2 = Entry()
        self.field_S2.pack()
        self.field_S3 = Entry()
        self.field_S3.pack()
        self.align_button = Button(master, text="Load from text fields", command=self.load_fields)
        self.align_button.pack()

        self.align_button = Button(master, text="DCA", command=self.align_dca_edu)
        self.align_button.pack()
        self.align_button = Button(master, text="Star", command=self.align_star)
        self.align_button.pack()
        self.align_button = Button(master, text="Progressive NJ", command=self.align_progressive_nj)
        self.align_button.pack()

        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()
        self.sequences = []
        for seq_record in SeqIO.parse("input.fasta", "fasta"):
            self.sequences.append((seq_record.description, seq_record.seq))

    def load_fasta(self, filename=""):
        """Loads sequences to be aligned from a fasta file"""

        if filename == "":
            self.sequences = []
            for seq_record in SeqIO.parse(self.field_input_file_name.get() + ".fasta", "fasta"):
                self.sequences.append((seq_record.description, seq_record.seq))
            print("Sequences loaded from file:")
            print(self.sequences)
        else:
            self.sequences = []
            for seq_record in SeqIO.parse(filename, "fasta"):
                self.sequences.append((seq_record.description, seq_record.seq))
            print("Sequences loaded from file:")
            print(self.sequences)
        return

    def load_fields(self):
        """Loads sequences to be aligned from the text fields"""
        self.sequences = []
        self.sequences.append(("Seq1", self.field_S1.get()))
        self.sequences.append(("Seq2", self.field_S2.get()))
        self.sequences.append(("Seq3", self.field_S3.get()))
        print("Sequences loaded from fields:")
        print(self.sequences)

    def align_dca(self):
        return



    def align_dca_edu(self, match_score = 1, mismatch_penalty = -1, gap_penalty = -1, extension_penalty = -1):
        """Performs DCA for 2 sequences"""
        msa = dca(list([self.sequences[0][1], self.sequences[1][1]]), 30, match_score, mismatch_penalty, gap_penalty, extension_penalty)
        print("ALIGNMENT:", "\n")
        for x in msa:
            print(x)
        save_msa_to_file(msa)
        return



    def align_star(self, match_score = 1, mismatch_penalty = -1, gap_penalty = -1, extension_penalty = -1, filename = "output.txt"):
        """Performs multiple sequence alignment by the center star method"""
        print(self.sequences)
        def extend(msa_to_extend, central, aligned_seq):
            """given the sequence aligned to the center sequence, adds the aligned sequence to the msa"""
            symbols_count = len(central) - central.count("-")
            sequence_to_append = ""
            msa_pointer = 0
            central_pointer = 0
            for i in range(symbols_count + 1):
                # iterate for each string of gaps, possibly empty
                # gap counters for current central sequence of the MSA and the central sequence of the pairwise alignment
                msa_counter = 0
                central_counter = 0
                if msa_pointer < len(msa_to_extend[0]):
                    while msa_to_extend[0][msa_pointer + msa_counter] == "-":
                        msa_counter += 1
                        if msa_pointer + msa_counter == len(msa_to_extend[0]):
                            break

                if central_pointer < len(central):
                    while central[central_pointer + central_counter] == "-":
                        central_counter += 1
                        if central_pointer + central_counter == len(central):
                            break


                if msa_counter == central_counter:
                    # same amount of gaps
                    sequence_to_append += aligned_seq[central_pointer:central_pointer + central_counter + 1]
                    central_pointer = central_pointer + central_counter + 1
                    msa_pointer = msa_pointer + msa_counter + 1
                elif msa_counter > central_counter:
                    # introduce gaps into the sequence being added to the msa
                    diff = msa_counter - central_counter
                    sequence_to_append += "-"*diff + aligned_seq[central_pointer:central_counter+central_pointer+1]
                    msa_pointer = msa_pointer + msa_counter + 1
                    central_pointer = central_pointer + central_counter + 1
                else:
                    # introduce gaps into the msa
                    diff = central_counter - msa_counter
                    for i, sequence in enumerate(msa_to_extend):
                        msa_to_extend[i] = msa_to_extend[i][:msa_pointer] + "-" * diff + msa_to_extend[i][msa_pointer:]
                    sequence_to_append += aligned_seq[central_pointer:central_pointer+central_counter+1]
                    central_pointer = central_pointer + central_counter + 1
                    msa_pointer = msa_pointer + msa_counter + diff + 1

            msa_to_extend.append(sequence_to_append)

            return

        #find the central sequence by pairwise alignments
        print("Finding the central sequence...")
        matrix = np.zeros((len(self.sequences), len(self.sequences)+1))
        for i, row in enumerate(matrix):
            for j in range(len(self.sequences)):
                if i == j:
                    continue
                matrix[i][j] = pairwise2.align.globalms(self.sequences[i][1], self.sequences[j][1],
                                                        match_score, mismatch_penalty, gap_penalty,
                                                        extension_penalty, score_only = True)
            matrix[i][len(self.sequences)] = np.sum(matrix[i][0:len(self.sequences)])
        central_sequence_score = matrix[0][-1]
        central_sequence = 0
        for i, row in enumerate(matrix):
            if row[len(row)-1] > central_sequence_score:
                central_sequence = i
                central_sequence_score = row[len(row)-1]
        print("CENTRAL SEQUENCE", central_sequence)
        print(str(self.sequences[central_sequence][1]))

        #obtain pairwise alignments with the central sequence
        print("Performing pairwise alignments with the central sequence...")
        alignments = []
        for i, sequence in enumerate(self.sequences):
            if i == central_sequence:
                print("Central sequence, skip", i, central_sequence)
                continue
            alignments.append(pairwise2.align.globalms(self.sequences[i][1],
                                                       self.sequences[central_sequence][1],
                                                        match_score, mismatch_penalty, gap_penalty,
                                                        extension_penalty, one_alignment_only = True)[0])
        print("Number of alignments obtained: ", len(alignments))
        print(alignments)
        # construct MSA iteratively by extending it with
        # every sequence other than the central sequence
        print("Creating the MSA...")
        msa = [str(self.sequences[central_sequence][1])]
        for alignment in alignments:
            central_pattern = alignment[1]
            insertion_pattern = alignment[0]
            extend(msa, central_pattern, insertion_pattern)
        #print the result to the console
        print("ALIGNMENT:")
        for x in msa:
            print(x)
        score = save_msa_to_file(msa, filename)
        return score

    def align_progressive_nj(self, match_score = 1, mismatch_penalty = -1, gap_penalty = -1, extension_penalty = -1, filename = "output.txt"):
        nodes_list = [Node([str(seq[1])]) for seq in self.sequences]

        calculator = DistanceCalculator('blosum62')

        distance_matrix = np.zeros((len(self.sequences), len(self.sequences)))

        for c in combinations(range(len(nodes_list)), 2):
            alignment = pairwise2.align.globalms(nodes_list[c[0]].consensus,
                                                 nodes_list[c[1]].consensus,
                                                 match_score, mismatch_penalty, gap_penalty,
                                                 extension_penalty, one_alignment_only=True)[0]

            aln = MultipleSeqAlignment([SeqIO.SeqRecord(Seq(alignment[0], Gapped(IUPAC.extended_protein, "-")), id="0"),
                                        SeqIO.SeqRecord(Seq(alignment[1], Gapped(IUPAC.extended_protein, "-")),
                                                        id="1")],
                                       Gapped(IUPAC.extended_protein, "-"))
            dm = calculator.get_distance(aln)
            distance_matrix[c[0]][c[1]] = distance_matrix[c[1]][c[0]] = dm[0][1]
        argmin = (0, 1)
        minvalue = distance_matrix[argmin[0], argmin[1]]
        for c in combinations(range(len(nodes_list)), 2):
            if distance_matrix[c[0]][c[1]] < minvalue:
                minvalue = distance_matrix[c[0]][c[1]]
                argmin = c
        print("ARGMIN, MIN", argmin, distance_matrix[argmin[0]][argmin[1]])

        print(distance_matrix)
        while len(nodes_list) > 1:
            argmin = (0, 1)
            minvalue = distance_matrix[argmin[0], argmin[1]]
            for c in combinations(range(len(nodes_list)), 2):
                if distance_matrix[c[0]][c[1]] < minvalue:
                    minvalue = distance_matrix[c[0]][c[1]]
                    argmin = c
            first = argmin[0]
            second = argmin[1]
            newnode = merge_nodes(nodes_list[first], nodes_list[second])
            nodes_list = nodes_list[0:first] + nodes_list[first + 1:second] + nodes_list[second + 1:]
            nodes_list.append(newnode)

            distance_matrix = np.zeros((len(nodes_list), len(nodes_list)))

            for c in combinations(range(len(nodes_list)), 2):
                alignment = pairwise2.align.globalms(nodes_list[c[0]].consensus,
                                                     nodes_list[c[1]].consensus,
                                                     match_score, mismatch_penalty, gap_penalty,
                                                     extension_penalty, one_alignment_only=True)[0]

                aln = MultipleSeqAlignment(
                    [SeqIO.SeqRecord(Seq(alignment[0], Gapped(IUPAC.extended_protein, "-")), id="0"),
                     SeqIO.SeqRecord(Seq(alignment[1], Gapped(IUPAC.extended_protein, "-")), id="1")],
                    Gapped(IUPAC.extended_protein, "-"))
                dm = calculator.get_distance(aln)
                distance_matrix[c[0]][c[1]] = distance_matrix[c[1]][c[0]] = dm[0][1]


        print("ALIGNMENT:")
        for x in nodes_list[0].msa:
            print(str(x))
        score = save_msa_to_file(nodes_list[0].msa, filename)
        return score


root = Tk()
gui = MSA(root)
root.mainloop()
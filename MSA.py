from tkinter import Tk, Label, Button, Entry, Checkbutton, IntVar
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo as mi
import numpy as np
import tests
class MSA:
    #sequences is a list of tuples (sequence_description, sequence)
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

    def load_fasta(self):
        """Loads sequences to be aligned from a fasta file"""

        self.sequences = []
        for seq_record in SeqIO.parse(self.field_input_file_name.get() +".fasta", "fasta"):
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
        return

    def align_dca(self):
        return

    def align_dca_edu(self):
        #dca for 2 sequences
        return



    def align_star(self, match_score = 1, mismatch_penalty = -1, gap_penalty = -1, extension_penalty = -1):
        """performs multiple sequence alignment by the center star method"""
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
                    # introduce gaps into msa

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
        msa = [self.sequences[central_sequence][1]]
        for alignment in alignments:
            central_pattern = alignment[1]
            insertion_pattern = alignment[0]
            extend(msa, central_pattern, insertion_pattern)
        #print the result to the console
        print("ALIGNMENT:")
        for x in msa:
            print(x)
        msa_array = np.array(msa)
        #save the result to file
        score = tests.compute_sum_of_pairs(msa_array, mi.blosum62)
        print(score)
        file = open("output.txt", "w")
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

    def align_progressive_nj(self):
        return



root = Tk()
gui = MSA(root)
root.mainloop()
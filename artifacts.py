#  this file contains code that may be useful,
#  but is not part of the final version of the project
from Bio.Align.Applications import TCoffeeCommandline
import os
from Bio.Align.Applications import ClustalOmegaCommandline
""" member function of the class MSA, used to evaluate implemented algorithms using given testcases,
    as well as generate commandline scripts to run clustal and tcoffee"""
# def compare_methods(self):
#     # filenames = ["testcase3_10",
#     #              "testcase3_20",
#     #              "testcase3_30",
#     #              "testcase3_40",
#     #              "testcase3_50",
#     #              "testcase3_100",
#     #              "testcase3_250",
#     #              "testcase3_500",
#     #              "testcase3_750",
#     #              "testcase3_1000",
#     #              "testcase10_10",
#     #              "testcase10_100",
#     #              "testcase10_250",
#     #              "testcase10_500",
#     #              "testcase10_750",
#     #              "testcase10_1000",
#     #              "testcase30_10",
#     #              "testcase30_100",
#     #              "testcase30_250",]
#     filenames = ["testcase30_2000",
#                  "testcase50_10",
#                  "testcase50_100",
#                  "testcase50_250",
#                  "testcase50_500"]
#     # "testcase50_1000",
#     # "testcase100_10",
#     # "testcase100_100",
#     # "testcase100_250",
#     #  ]
#     with open("SCORES", 'w') as file:
#         file2 = open("run_testcases_tcoffee", "w")
#         file3 = open("run_testcases_clustalo", "w")
#         for filename in filenames:
#             print("COMPUTING ", filename)
#             self.load_fasta(filename + ".fasta")
#             star_score = self.align_star(filename=filename + "_output_star")
#             progressive_nj_score = self.align_progressive_nj(filename=filename + "_output_progressive_nj")
#             file.write(filename)
#             file.write("\n")
#             file.write("STAR SCORE: ")
#             file.write(str(star_score))
#             file.write(", PROGRESSIVE NJ SCORE: ")
#             file.write(str(progressive_nj_score))
#             file.write("\n")
#
#             command = str(TCoffeeCommandline(infile=filename,
#                                              output="clustalw",
#                                              outfile=filename + "_output_tcoffee" + ".aln"))
#             file2.write(command)
#             file2.write("\n")
#             os.system(command)
#
#             command = str(ClustalOmegaCommandline(infile=filename + ".fasta",
#                                                   outfile=filename + "_output_clustalo" + ".aln"))
#             file3.write(command)
#             file3.write("\n")
#             os.system(command)
#         file2.close()
#         file3.close()

from Bio.Align.Applications import TCoffeeCommandline
filenames = ["testcase3_10",
                     "testcase3_20",
                     "testcase3_30",
                     "testcase3_40",
                     "testcase3_50",
                     "testcase3_100",
                     "testcase3_250",
                     "testcase3_500",
                     "testcase3_750",
                     "testcase3_1000",
                     "testcase10_10",
                     "testcase10_100",
                     "testcase10_250",
                     "testcase10_500",
                     "testcase10_750",
                     "testcase10_1000",
                     "testcase30_10",
                     "testcase30_100",
                     "testcase30_250",]

# file2 = open("run_testcases_tcoffee", "w")
# for filename in filenames:
#     command = str(TCoffeeCommandline(infile=filename+".fasta",
#                                      output="clustalw",
#                                      outfile=filename + "_output_tcoffee" + ".aln"))
#     file2.write(command)
#     file2.write("\n")
#
# file2.close()
#
# file2 = open("results_tcoffee", "w")
#
# for filename in filenames:
#     score = evaluate_clustal(filename + "_output_tcoffee.aln", compute_sum_of_pairs)
#     file2.write(filename)
#     file2.write("\n")
#     file2.write(str(score))
#     file2.write("\n")
#
# file2.close()

from Bio.Align.Applications import ClustalOmegaCommandline


file2 = open("run_testcases_clustalo", "w")
for filename in filenames:
    command = str(ClustalOmegaCommandline(infile=filename+".fasta",
                                     outfile=filename + "_output_clustalo" + ".aln"))
    file2.write(command)
    file2.write("\n")

file2.close()

# file2 = open("results_clustalo", "w")
#
# for filename in filenames:
#     score = evaluate_clustal(filename + "_output_clustalo.aln", compute_sum_of_pairs)
#     file2.write(filename)
#     file2.write("\n")
#     file2.write(str(score))
#     file2.write("\n")
#
# file2.close()

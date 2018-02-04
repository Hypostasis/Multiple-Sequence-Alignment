from tkinter import Tk, Label, Button, Entry, Checkbutton, IntVar
from Bio import SeqIO
class MSA:
    sequences = []
    def __init__(self, master):
        self.master = master
        master.title("Multiple Sequence Alignment")

        self.label = Label(master, text = "Enter your sequences here or as a FASTA file:")
        self.label.pack()
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


        self.align_button = Button(master, text="DCA", command=self.align_dca)
        self.align_button.pack()
        self.align_button = Button(master, text="Star", command=self.align_star)
        self.align_button.pack()
        self.align_button = Button(master, text="Progressive NJ", command=self.align_progressive_nj)
        self.align_button.pack()

        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()

        self.data = []

    def load_fasta(self):

        for seq_record in SeqIO.parse("input.fasta", "fasta"):
            self.sequences.append(seq_record.seq)
        return
    def load_fields(self):
        self.sequences.append(self.field_S1.get())
        self.sequences.append(self.field_S2.get())
        self.sequences.append(self.field_S3.get())
        return

    def align_dca(self):
        return

    def align_star(self):
        for seq in self.sequences:
            print(seq)
        return
    def align_progressive_nj(self):
        return



root = Tk()
gui = MSA(root)
root.mainloop()
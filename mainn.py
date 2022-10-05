import operator
import random
import streamlit as st

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt

from Suffix_Tree import SuffixTree


class Biopython:

    # SEARCH_PATTERNS

    def Search_Patterns(self):
        self.full_seq = input("Enter your seq: ")
        patterns = list(
            map(str, input("Enter the patterns separated by a space: ").split()))

        Trie = SuffixTree(self.full_seq)
        for i in patterns:
            if Trie.hasSubstring(i):
                print(f"pattern \"{i}\": found")
            else:
                print(f"pattern \"{i}\": not found")

# SEQ_GENERATOR

    def new_rnd_seq(sl):
        s = ''
        for x in range(sl):
            s += random.choice('ATCG')
        # s += random.sample(’ATCG’,1)[0] is not so fast.
        return s

    def Generate_fasta(self, seqlen1, seqlen2, seq_num):
        newfh = open("randomseqs.txt", "w+")

        for i in range(1, seq_num + 1):
            rsl = random.randint(min(seqlen1, seqlen2), max(seqlen1, seqlen2))

            # Generate the random sequence
            rawseq = Biopython.new_rnd_seq(rsl)

            # Generate a correlative name
            seqname = "Sequence_number_" + str(i)

            rec = SeqRecord(Seq(rawseq), id=seqname, description="")
            SeqIO.write([rec], newfh, "fasta")

        newfh.close()

        f = open('randomseqs.txt', 'r')
        file_contents = f.read()
        return file_contents


# =====================================================================================
# FREQUENT_PATTERNS

    def Frequent_Patterns(self, seq, K_mer):
        self.seq = seq

        Frequent_Patterns = {}
        counts = []
        for i in range(0, len(self.seq) - K_mer + 1):
            pattern = self.seq[i: (i + K_mer)]
            if pattern not in Frequent_Patterns:
                Frequent_Patterns[pattern] = self.seq.count(pattern)

        Max = max(Frequent_Patterns.items(), key=operator.itemgetter(1))[1]
        return ", ".join([k for k, v in Frequent_Patterns.items() if v == Max])


# =====================================================================================
# MINIMUM_SKEW

    def count_GC(self, seq):
        sequance = seq
        C = 0
        G = 0

        for i in sequance:
            if i == 'G':
                G += 1

            elif i == 'C':
                C += 1

        length = len(sequance)
        if length != 0:
            GC = (G / length) + (C / length)*100
            roundGC = str(round(GC, 2))
            fulltxt = "The GC content: " + roundGC + " %"
            return fulltxt

    def Minimum_Skew(self, seq):

        self.full_seq = seq
        SKEW_list = [0]
        SKEW = 0
        for i in self.full_seq:
            if i == 'C':
                SKEW -= 1
            if i == "G":
                SKEW += 1
            SKEW_list.append(SKEW)

        # The skew diagram for Genome
        y = SKEW_list
        x_labels = [i for i in self.full_seq]
        x_labels.insert(0, 0)

        plt.plot(range(len(x_labels)), y)
        plt.xticks(range(len(x_labels)), x_labels)

        plt.show()
################################################

    def transcript(x):

        l = str(x)
        l = l.replace("C", "G")
        l = l.replace("G", "C")
        l = l.replace("A", "U")
        l = l.replace("T", "A")
        return l

    if __name__ == "__main__":
        x = input("enter seq  ")
        print(transcript(x.upper()))

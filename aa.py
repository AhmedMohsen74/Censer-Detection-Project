from typing import Pattern
import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image
import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
from mainn import Biopython
from breastCancer import model
from Suffix_Tree import SuffixTree

nav = st.sidebar.radio("Choose operation",
                       ["Breast cancer detection", "Search for patterns", "Find frequent pattern", "Generate", "Minimum Skew", "Transcription"])
#################################
if nav == "Breast cancer detection":

    side, beside = st.columns(2)
    with side:
        st.title("Detecting Breast Cancer")
    with beside:
        st.image("detect.png")

    seq1 = st.number_input("Clump Thikness: ", key="seq1")
    seq2 = st.number_input("Uniform Cell Size: ", key="seq2")
    seq3 = st.number_input("Uniform Cell Shape: ", key="seq3")
    seq4 = st.number_input("Marginal Adhesion: ", key="seq4")
    seq5 = st.number_input("Single Epithelial Cell Size: ", key="seq5")
    seq6 = st.number_input("Bare Nuclei: ", key="seq6")
    seq7 = st.number_input("Bland Chromatin: ", key="seq7")
    seq8 = st.number_input("Normal Nucleoli: ", key="seq8")
    seq9 = st.number_input("Mitosis: ", key="seq9")

    if(st.button("Predict Cancer")):
        input_features = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9]
        features_value = [np.array(input_features)]
        features_name = ['clump_thickness', 'uniform_cell_size', 'uniform_cell_shape',
                         'marginal_adhesion', 'single_epithelial_size', 'bare_nuclei',
                         'bland_chromatin', 'normal_nucleoli', 'mitoses']
        df = pd.DataFrame(features_value, columns=features_name)
        output = model.predict(df)
        if output == 4:
            st.write("Patient Has Breast cancer")
        else:
            st.write("Patient Has No Breast cancer")


#################################
if nav == "Search for patterns":
    side, beside = st.columns(2)
    with side:
        st.title('Search for patterns')
    with beside:
        st.image("SearchForPatterns.png")

    seq = st.text_input('Enter your seq',
                        '').upper()

    patt = st.text_input('Enter your pattern',
                         '').upper()

    patterns = patterns = list(map(str, patt.split()))

    Trie = SuffixTree(seq)
    if(st.button("Start The Process")):
        for i in patterns:
            if Trie.hasSubstring(i):
                st.write("Pattern Found")
            else:
                st.write("Pattern Not Found")
#################################
if nav == "Find frequent pattern":

    side, beside = st.columns(2)
    with side:
        st.title('Find Frequent Patterns')
    with beside:
        st.image('find.png')

    seq = st.text_input("Enter your sequence: ", key="seq")
    k_mer = st.number_input("Enter the k-mer: ", key="k_mer", step=1)
    freq_patterns = Biopython.Frequent_Patterns(
        self=Biopython, seq=seq, K_mer=k_mer)
    if(st.button("Start The Process")):
        st.write(freq_patterns)

#################################

if nav == "Generate":

    newfh = open('randomseqs.fasta', 'w')
    seqNum = st.number_input("Enter Number of Seq", 1, 1000, 500, 1)
    st.text("Numb of nitrogen bases")

    f = st.number_input("From", 1, 10000, 4000, 1)  # from
    t = st.number_input("To", 1, 20000, 15000, 1)  # to

    if st.button("Generate"):
        if(f < 80 and t >= 200):
            st.write("You Must Enter in range from 4000 , 15000")
        else:
            for i in range(1, seqNum):
                # Creates a random number in the range of 4000-15000
                rsl = random.randint(f, t)
                rseq = Biopython.new_rnd_seq(rsl)  # raw sequance
                seqn = 'SequenceNumber' + str(i)  # seqance name
                rec = SeqRecord(Seq(rseq), id=seqn, description='')
                SeqIO.write([rec], newfh, 'fasta')

            with open("randomseqs.fasta") as ran:
                st.download_button("Download the Sequence file",
                                   data=ran, file_name=("randomseqs.fasta"))


#################################
if nav == "Minimum Skew":

    side, beside = st.columns(2)
    with side:
        st.title('Minimun Skew')
    with beside:
        st.image('gene.jpeg')

    seq = st.text_input("Enter your sequence: ", key="seq")
    freq_patterns = Biopython.Minimum_Skew(self=Biopython, seq=seq)
    if(st.button("Start The Process")):
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot(freq_patterns)
        st.write(Biopython.count_GC(self=Biopython, seq=seq))

#################################

if nav == "Transcription":
    side, beside = st.columns(2)
    with side:
        st.title('Enter Your Pattern')
    with beside:
        st.image("SearchForPatterns.png")

    seq = st.text_input('Enter your seq', '').upper()

    if(st.button("Start The Process")):
        Trie = Biopython.transcript(seq)
        st.write(Trie)

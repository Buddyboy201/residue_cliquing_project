import atom
import residue
import protein
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt

class ProteinDatasetAnalyzer:
    def __init__(self, proteins=[]):
        self.proteins = proteins

    def add_protein(self, new_protein):
        self.proteins.append(new_protein)

    def get_proteins(self):
        return self.proteins

    def get_protein(self, name):
        x = protein.binary_search(self.proteins, name)
        if x == -1: return -1
        return self.proteins[x]

    #protein bench data format:
    # [ [name (same as in private field), clique file path, index file path], ...]

    def get_general_clique_stats(self, bench_protein_data=None):
        for P in self.proteins:
            print("{} GENERAL CLIQUE STATS".format(P.get_name()))
            P.get_cliques("centroid")
            P.get_cliques("atom")
            protein.distance_analysis("atom")
            if bench_protein_data != None:
                for N, C, I in bench_protein_data:
                x = protein.binary_search(self.proteins
            protein.freq_analysis("atom", test_bench_cliques_file)
            protein.bench_to_new_analysis("atom", test_bench_cliques_file, test_bench_index_file)

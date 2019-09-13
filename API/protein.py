import atom
import residue
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt



#generate clique codes
#sort clique codes
#find included bench clique codes in new cliques of any type

def binary_search(arr, val):
    start = 0
    end = len(arr)-1
    while True:
        mid = int((start+end)/2)
        #print(start, end)
        if start >= end and arr[mid] != val: return -1
        if arr[mid] == val: return mid
        elif val > arr[mid]: start = mid+1
        else: end = mid-1

###############################################################################################
#TODOs:
#   make option to exclude/include backbone atoms(50% completion)
#   csv data logging
#   matplotlib/excel/google_sheets graphing methods
#   other data analysis methods[*priority*: benchmark to new analysis]
###############################################################################################

def get_dist(coord1, coord2):
    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)

class Protein:
    def __init__(self, name, file_path):
        self.name = name
        self.file_path = file_path
        self.residues = {}
        self.atom_cliques = None
        self.centroid_cliques = None
        #self.distance_cutoff = 6
        self.centroid_clique_frequency = None
        self.atom_clique_frequency = None
        with open(self.file_path) as pdb_file:
            for line in pdb_file:
                if line[0:4] == "ATOM":
                    res_id = int(line[22:26].strip(" "))
                    atom_id = int(line[6:11].strip(" "))
                    res_name = line[17:20].strip(" ")
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    symbol = line[76:78].strip(" ")
                    atom_name = line[12:16].strip(" ")
                    coords = (coordx, coordy, coordz)
                    atm = atom.Atom(symbol, atom_name, atom_id, coords)
                    if self.residues.get(res_id) == None:
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm])
                    else: self.residues[res_id].add_atom(atm)
        for i in self.residues:
            self.residues[i].update_COM()

    def convert_to_code(self, clique):
        clique = list(clique)
        clique.sort()
        s = ""
        l = list(map(str, clique))
        for i in l: s+= i
        return int(s)

    def get_sorted_clique_codes(self, cliques):
        codes = []
        for i in range(len(cliques)):
            #print(cliques[i])
            codes.append(self.convert_to_code(cliques[i]))
        codes.sort()
        return codes

    def get_resindex_resid_hash(self, bench_index_file):
        resindex_resid = {}
        with open(bench_index_file) as file:
            for line in file:
                a = line.strip(" \n").split(" ")
                resindex_resid[int(a[0])] = int(a[1])
        return resindex_resid

    def convert_resindex_resid(self, bench_cliques, resindex_resid):
        for i in range(len(bench_cliques)):
            for j in range(len(bench_cliques[i])):
                bench_cliques[i][j] = resindex_resid[bench_cliques[i][j]]
        return bench_cliques

    def get_bench_cliques(self, bench_file):
        cliques = []
        with open(bench_file) as file:
            for line in file:
                pos = line.find("|", line.find("|")+1, len(line)-1)
                clique = list(map(int, line[pos+2:].strip(" \n").split(" ")))
                cliques.append(clique)
        return cliques

    def get_included_bench_cliques(self, codes_bench, codes_atom):
        included = []
        for i in codes_bench:
            x = binary_search(codes_atom, i)
            if x == -1: continue #print("\n\n\n\n\n", x, "\n\n\n\n\n")
            else:
                included.append(codes_atom[x])
                #print("\n\n\n\n\n", x, "\n\n\n\n\n")
        return included
    
    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def get_cliques(self, clique_type, exclude_backbone=False, distance_cutoff=6):
        if clique_type == "centroid":
            if self.centroid_cliques == None:
                self.generate_cliques("centroid", exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
            return self.centroid_cliques
        elif clique_type == "atom":
            if self.atom_cliques == None:
                self.generate_cliques("atom", exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
            return self.atom_cliques
        else: raise Exception("invalid clique type")
            
    def generate_centroid_cliques(self, exclude_backbone=False, distance_cutoff=6):
        centroids = []
        centroid_res = {}
        for i in self.residues:
            centroids.append(self.residues[i].get_centroid())
            centroid_res[self.residues[i].get_centroid()] = self.residues[i]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)

    def generate_atom_cliques(self, exclude_backbone=False, distance_cutoff=6):
        coords = []
        coords_resid = {}
        for i in self.residues:
            for j in self.residues[i].get_atoms():
                coords.append(j.get_coords())
                coords_resid[j.get_coords()] = i
        coords_array = np.array(coords)
        del coords
        tri = scipy.spatial.qhull.Delaunay(coords_array)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.atom_cliques = list(nx.find_cliques(graph))
        temp_arr = []
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = coords_resid[tuple(list(coords_array[self.atom_cliques[i][j]]))]
            self.atom_cliques[i] = list(set(self.atom_cliques[i]))
            self.atom_cliques[i].sort()
            if self.atom_cliques[i] not in temp_arr:
                temp_arr.append(self.atom_cliques[i])
        self.atom_cliques = np.array(temp_arr)
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = self.residues[self.atom_cliques[i][j]]

    def generate_cliques(self, clique_type, exclude_backbone=False, distance_cutoff=6):
        if clique_type == "centroid":
            self.generate_centroid_cliques(exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
        elif clique_type == "atom":
            self.generate_atom_cliques(exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff)
        else: raise Exception("Invalid clique type")

    def getMaxMinDistance(self, coords):
        max_dist = 0
        min_dist = 10000
        if len(coords) == 1: return 0, 0
        for i in range(len(coords)-1):
            for j in range(i+1, len(coords)):
                dist = get_dist(coords[i], coords[j])
                max_dist = max([max_dist, dist])
                min_dist = min([min_dist, dist])
        return max_dist, min_dist

    def distance_analysis(self, clique_type):
        if clique_type == "centroid":
            print("Min/Max CENTROID CLIQUE STATS")
            clique_max_sum = 0
            clique_min_sum = 0
            count = 0
            for i in self.centroid_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                if clique_max != 0 and clique_min != 0:
                    clique_max_sum += clique_max
                    clique_min_sum += clique_min
                    count += 1
            clique_max_avg = float(clique_max_sum)/count
            clique_min_avg = float(clique_min_sum)/count
                #print("min_dist: {} | max_dist: {}".format(clique_min, clique_max))
            print("avg_min_dist: {} | avg_max_dist: {}".format(clique_min_avg, clique_max_avg))
        elif clique_type == "atom":
            print("Min/Max ATOM CLIQUE STATS")
            clique_max_sum = 0
            clique_min_sum = 0
            count = 0
            for i in self.atom_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                if clique_max != 0 and clique_min != 0:
                    clique_max_sum += clique_max
                    clique_min_sum += clique_min
                    count += 1
            clique_max_avg = float(clique_max_sum)/count
            clique_min_avg = float(clique_min_sum)/count
                #print("min_dist: {} | max_dist: {}".format(clique_min, clique_max))
            print("avg_min_dist: {} | avg_max_dist: {}".format(clique_min_avg, clique_max_avg))
        else: raise Exception("Invalid clique type")

    def freq_analysis(self, clique_type, bench_cliques_file=None):
        freq_arr = [0,0,0,0,0,0,0]
        if clique_type == "centroid":
            for i in self.centroid_cliques:
                freq_arr[len(i)] += 1
            self.centroid_clique_frequency = freq_arr
        elif clique_type == "atom":
            for i in self.atom_cliques:
                freq_arr[len(i)] += 1
            self.atom_clique_frequency = freq_arr
        else: raise Exception("Invalid clique type")
        bench_freq_arr = [0,0,0,0,0,0,0]
        if bench_cliques_file != None:
            bench_cliques = self.get_bench_cliques(bench_cliques_file)
            for i in bench_cliques:
                bench_freq_arr[len(i)] += 1
            print("{} v. BENCH CLIQUE SIZE FREQS".format(clique_type))
            print("{}: {}".format(clique_type, freq_arr[1:]))
            print("bench: {}".format(bench_freq_arr[1:]))
        else:
            print("{} CLIQUE SIZE FREQS".format(clique_type))
            print("{}: {}".format(clique_type, freq_arr[1:]))
        #return freq_arr

    def get_included_bench_cliques(self, codes_bench, codes_new):
        included = []
        for i in codes_bench:
            x = binary_search(codes_new, i)
            if x != -1:
                included.append(codes_new[x])
        return included

    def bench_to_new_analysis(self, clique_type, bench_cliques_file, bench_index_file):
        resindex_resid = self.get_resindex_resid_hash(bench_index_file)
        cliques = []
        new_codes = None
        if clique_type == "centroid":
            for i in self.centroid_cliques:
                clique = []
                for j in i:
                    clique.append(j.get_resid())
                cliques.append(clique)
            cliques = np.array(cliques)
            new_codes = self.get_sorted_clique_codes(cliques)
        elif clique_type == "atom":
            for i in self.atom_cliques:
                clique = []
                for j in i:
                    clique.append(j.get_resid())
                cliques.append(clique)
            cliques = np.array(cliques)
            new_codes = self.get_sorted_clique_codes(cliques)
        else: raise Exception("Invalid clique type")
        bench_cliques = self.get_bench_cliques(bench_cliques_file)
        bench_cliques = self.convert_resindex_resid(bench_cliques, resindex_resid)
        bench_codes = self.get_sorted_clique_codes(bench_cliques)
        included = self.get_included_bench_cliques(bench_codes, new_codes)
        #print(included)
        print("{} CLIQUES BENCH INCLUSION STATS".format(clique_type)) 
        print("% included: {}".format((float(len(included))/len(bench_codes))*100))
        print("New: {} | Bench: {} | Included: {}".format(len(new_codes), len(bench_codes), len(included)))
        
test_bench_cliques_file = "C:\\pdb_data\\domains\\f1a4fa_.nomc.cliques"
test_bench_index_file = "C:\\pdb_data\\domains\\f1a4fa_.index"

protein = Protein("f1a4fa_", "C:\\pdb_data\\domains\\f1a4fa_.pdb")
protein.get_cliques("centroid")
protein.get_cliques("atom")
protein.distance_analysis("atom")
protein.freq_analysis("atom", test_bench_cliques_file)
protein.bench_to_new_analysis("atom", test_bench_cliques_file, test_bench_index_file)

#print(protein.centroid_cliques[0][1].get_centroid())


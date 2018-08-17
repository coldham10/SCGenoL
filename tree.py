import numpy as np

class Tree:
    """A tree of cells constructed by neighbour joining
    by Hamming distance."""

    def __init__(self, loci, ado_rate): #fp_rate, cyto_deam_rate):
        self.prior_ado_rate = ado_rate
        self.loci           = loci
        self.n_cells        = len(loci[0].n_cells)
        self.dist_matrix    = np.zeros(self.n_cells, self.n_cells)
        self.hamming_dist()
        self.tree()

    def hamming_dist(self):
        #TODO
        # Is hamming distance correct for a probabilistic tree? weighted by likelihoods?
        pass

    def tree(self):
        #TODO
        pass

    def estimate_ado_rate(self):
        #TODO
        pass

    #estimate fp_rate, cyto_deam_rate, CNV_rate, SFS, 

class Node:
    """A node on the tree consisting of a single cell"""
    def __init__(self, cell_no, parent):
        self.cell_no    = cell_no
        self.parent     = parent
        self.children   = []


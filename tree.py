import numpy as np

class Tree:
    """A tree of cells constructed by neighbour joining
    by Hamming distance."""

    def __init__(self, loci, ado_rate): #fp_rate, cyto_deam_rate):
        self.prior_ado_rate = ado_rate
        self.loci           = loci
        self.n_cells        = loci[0].n_cells
        self.dist_matrix    = np.zeros((self.n_cells, self.n_cells))
        self.hamming_dist()
        self.tree()

    def hamming_dist(self):
        # Is hamming distance correct for a probabilistic tree? weighted by posterior probabilities?
        for locus in self.loci:
            for i in range(self.n_cells):
                for j in range(i, self.n_cells):
                    # infinite sites, considering only g=0 or g=1
                    p01 = locus.cells[i].log_probs[0] + locus.cells[j].log_probs[1]
                    p10 = locus.cells[i].log_probs[1] + locus.cells[j].log_probs[0]
                    locus_dist = np.logaddexp(p01, p10)
                    prev_dist = self.dist_matrix[i,j]
                    dist = np.logaddexp(prev_dist, locus_dist)
                    self.dist_matrix[i,j] = dist
                    self.dist_matrix[j,i] = dist

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


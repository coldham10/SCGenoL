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
        for locus in self.loci:
            for i in range(self.n_cells):
                for j in range(i, self.n_cells):
                    # infinite sites, considering only g=0 or g=1
                    p01 = locus.cells[i].log_probs[0] + locus.cells[j].log_probs[1]
                    p10 = locus.cells[i].log_probs[1] + locus.cells[j].log_probs[0]
                    #TODO: this assumes independence of the two cells. How is this justified?
                    locus_dist = np.logaddexp(p01, p10)
                    prev_dist = self.dist_matrix[i,j]
                    dist = np.logaddexp(prev_dist, locus_dist)
                    self.dist_matrix[i,j] = dist
                    self.dist_matrix[j,i] = dist

    def add_root_dists(self, dist_array):
        dist_array = np.insert(dist_array, 0, 0, axis=0)
        dist_array = np.insert(dist_array, 0, 0, axis=1)
        for locus in self.loci:
            for i in range(self.n_cells):
                p10 = locus.cells[i].log_probs[1] 
                prev_dist = dist_array[i,0]
                dist = np.logaddexp(prev_dist, p10)
                dist_array[i,0] = dist
                dist_array[0,i] = dist

    def calculate_q(self, d):
        n = d.shape[0]
        q = np.zeros((n,n))
        for i in range(n):
            for j in range(i+1,n):
                isum = np.sum(d[i, ...])
                jsum = np.sum(d[j, ...])
                qval = (n-2)*d[i,j] - isum - jsum
                q[i,j] = q[j,i] = qval
        return q

    def tree(self):
        #TODO test
        active_nodes = []
        distances = self.dist_matrix.copy()
        self.root = Node(False, cell_no=-1)
        active_nodes.append(root)
        self.add_root_dists(distances)
        for i in range(self.n_cells):
            node = Node(True, cell_no=i)
            active_nodes.append[node]
        while len(active_nodes) > 2:
            q_mat = self.calculate_q(distances)

    def estimate_ado_rate(self):
        #TODO
        pass

    #estimate fp_rate, cyto_deam_rate, CNV_rate, SFS, 

class Node:
    """A node on the tree consisting of a single cell"""
    def __init__(self, leaf, cell_no=None):
        self.neighbours = []
        self.is_leaf = leaf
        self.cell_no = cell_no

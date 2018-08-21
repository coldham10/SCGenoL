"""For building a tree based on posterior genotype probabilities.

Classes: Tree, Node

Tree(loci) is initialized with an ordered list of Locus objects. The 
probabilistic Hamming distances are computed (in log-space) for each
pair of cells, and the tree is subsequently built by neighbour joining.

The tree also includes a root node using the reference or germline (if 
provided) genotype. The tree is then rooted at this node."""

__author__ = "C. F. Oldham"

import numpy as np

class Tree:
    """A tree of cells constructed by neighbour joining by Hamming 
    distance."""

    root = None

    def __init__(self, loci) #ado_rate, fp_rate, cyto_deam_rate):
        self._loci = loci
        self._n_cells = loci[0].n_cells
        self._dist_matrix = np.zeros((self.n_cells, self.n_cells))
        self._init_hamming_dists()
        self._tree()

    def estimate_ado_rate(self):
        #TODO
        pass

    def _init_hamming_dists(self):
        """Initializes the log of probabilistic Hamming distances 
        between real cells."""
        for locus in self._loci:
            for i in range(self._n_cells):
                for j in range(i, self._n_cells):
                    # TODO: add g=2 to this metric
                    # infinite sites, considering only g=0 or g=1
                    p01 = locus.cells[i].log_probs[0] 
                          + locus.cells[j].log_probs[1]
                    p10 = locus.cells[i].log_probs[1] 
                          + locus.cells[j].log_probs[0]
                    #TODO: this assumes independence of cells. Justified?
                    locus_dist = np.logaddexp(p01, p10)
                    prev_dist = self._dist_matrix[i,j]
                    dist = np.logaddexp(prev_dist, locus_dist)
                    self._dist_matrix[i,j] = dist
                    self._dist_matrix[j,i] = dist

    def _tree(self):
        #TODO test
        distances = self._dist_matrix.copy()
        active_nodes = []
        self.root = Node(False, cell_no=-1)
        # Root treated as extra sample (pos 0) with heterozygous welltype
        active_nodes.append(root)
        self._add_root_dists(distances)
        for i in range(self._n_cells):
            node = Node(True, cell_no=i)
            active_nodes.append[node]
        while len(active_nodes) > 2:
            q_mat = self._calculate_q(distances)
            n1, n2 = np.unravel_index(q_mat.argmin(), q_mat.shape)
            new_node = Node(False)
            self._link_neighbours(new_node, active_nodes[n1], active_nodes[n2])
            d1 = distances[n1, ...]
            d2 = distances[n2, ...]
            self._deactivate(active_nodes, distances, [n1, n2])
            active_nodes.append(new_node)
            new_dists = self._new_node_dists(n1, n2, d1, d2)
            distances = self._update_dists(distances, new_dists)
        active_nodes[0].neighbours.append(active_nodes[1])
        active_nodes[1].neighbours.append(active_nodes[0])
        self._root_tree()


    def _add_root_dists(self, dist_array):
        """Adds a pseudo-cell of the root to a distance matrix at 
        index 0 """
        dist_array = np.insert(dist_array, 0, 0, axis=0)
        dist_array = np.insert(dist_array, 0, 0, axis=1)
        for locus in self._loci:
            for i in range(self._n_cells):
                p10 = locus.cells[i].log_probs[1] 
                prev_dist = dist_array[i,0]
                dist = np.logaddexp(prev_dist, p10)
                dist_array[i,0] = dist
                dist_array[0,i] = dist

    def _calculate_q(self, d):
        """Calculates the Q matrix used in NJ"""
        n = d.shape[0]
        q = np.zeros((n,n))
        for i in range(n):
            for j in range(i,n):
                if i == j:
                    # So argmin of neighbour joining doesn't join the same node
                    q[i,j] = np.inf
                else:
                    isum = np.sum(d[i, ...])
                    jsum = np.sum(d[j, ...])
                    qval = (n-2)*d[i,j] - isum - jsum
                    q[i,j] = q[j,i] = qval
        return q

    def _new_node_dists(self, i1, i2, d1, d2):
        """Returns a list of logged distances between each node and
        the new node created between i1 and i2 given distances from the
        old matrix"""
        d12 = d1[2]
        new_dists = []
        for k in range(len(d1)):
            if k != i1 and k != i2:
                new_dk = np.log(0.5) + np.logaddexp(d1[k], d2[k], -d12)
                new_dists.append(new_dk)
        return new_dists

    def _update_dists(self, dists, to_add):
        n = dists.shape[0]
        dists  = np.insert(dists, n, to_add, axis=0)
        to_add = np.insert(to_add, n, 0)
        dists  = np.insert(dists, n, to_add, axis=1)
        return dists

    def _link_neighbours(self, node, n1, n2):
        """Links the newly created 'node' with the two neighbours that
        joined to form it"""
        node.neighbours.append(n1)
        node.neighbours.append(n2)
        n1.neighbours.append(node)
        n2.neighbours.append(node)

    def _deactivate(self, active_list, dists, indices):
        """Removes nodes given in indices from the list of active 
        nodes, as well as pruning the distance matrix accordingly"""
        for node in indeces:
            del active_list[node]
            np.delete(dists, node, axis=0)
            np.delete(dists, node, axis=1)


    #estimate fp_rate, cyto_deam_rate, CNV_rate, SFS, 

class Node:
    """A node on the tree consisting of a single cell"""
    def __init__(self, is_leaf, cell_no=None):
        self.neighbours = []
        self.children   = []
        self.is_leaf    = leaf
        self.cell_no    = cell_no
        # For rooting
        self.visited    = False

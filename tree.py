class tree:
    """A tree of cells constructed by neighbour joining
    by Hamming distance."""

    def __init__(self, loci, ado_rate): #fp_rate, cyto_deam_rate):
        self.prior_ado_rate = ado_rate
        self.loci = loci
        l1 = loci[0]

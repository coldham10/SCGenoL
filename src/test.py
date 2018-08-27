import parser
import numpy as np
import tree
def test():
    pfile = parser.Pileup_file("test2.pileup")
    amp_mat = parser.Amplification_matrix(fp_error = 0.01)
    loci = []
    for locus in pfile:
        loci.append(locus)
        for c in locus.cells:
            c.calculate_naive_posteriors(amp_mat.matrix,0.01,0.0001)
    
    t = tree.Tree(loci)
    return (loci, t)

if __name__ == "__main__":
    l, t = test()
    t.show()

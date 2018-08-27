import parser
import numpy as np
import tree
import time
def test():
    pfile = parser.Pileup_file("tmp.pileup")
    amp_mat = parser.Amplification_matrix(fp_error = 0.01)
    loci = []
    locus_init_time = 0
    posterior_time  = 0 
    tree_time       = 0
    last_time = time.perf_counter()
    for locus in pfile:
        loci.append(locus)
        now = time.perf_counter()
        locus_init_time += now - last_time
        last_time = now
        for c in locus.cells:
            c.calculate_naive_posteriors(amp_mat.matrix,0.01,0.0001)
        now = time.perf_counter()
        posterior_time += now - last_time
        last_time = now

    t = tree.Tree(loci)
    now = time.perf_counter()
    tree_time += now - last_time
    return (locus_init_time, posterior_time, tree_time, t)

if __name__ == "__main__":
    t1, t2, t3, t = test()
    t.show()
    out = "Locus init: {0}, Posterior: {1}, tree: {2}".format(t1, t2, t3)
    print(out)

import parser
import numpy as np
pfile = parser.Pileup_file("test.pileup")
amp_mat = parser.Amplification_matrix(fp_error = 0.01)
for locus in pfile:
    probs = []
    for c in locus.cells:
        c.calculate_naive_posteriors(amp_mat.matrix,0.01,0.0001)
        probs.append(np.exp(c.log_probs))

    #print("{0:3.3f}\t{0:3.3f}\t{0:3.3f}\t".format(probs[0][0],probs[0][1],probs[0][2]))
    print(probs)

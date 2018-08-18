import parser
import numpy as np
pfile = parser.Pileup_file("test.pileup")
l = next(pfile)
l = next(pfile)
c = l.cells[0]
amp_mat = parser.Amplification_matrix(fp_error = 0.01)
c.calculate_naive_posteriors(amp_mat.matrix,0.01,0.0001)
print (np.exp(c.log_probs))

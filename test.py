import parser
import numpy as np
pfile = parser.Pileup_file("test.pileup")
l = next(pfile)
l = next(pfile)
c = l.first_cell
amp_mat = parser.Amplification_matrix(fp_error = 0.01)


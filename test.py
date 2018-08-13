import parser
import numpy as np
pfile = parser.Pileup_file("test.pileup")
l = next(pfile)
l = next(pfile)
c = l.first_cell
amp_mat = parser.Amplification_matrix(fp_error = 0.01)
print (np.exp(c.l_genotype_from_SC_reads(0,amp_mat.matrix,0.002)))
print (np.exp(c.l_genotype_from_SC_reads(1,amp_mat.matrix,0.002)))
print (np.exp(c.l_genotype_from_SC_reads(2,amp_mat.matrix,0.002)))

import sys
from locus import Locus
from collections import Iterator

class Pileup_file(Iterator):
    """Class implements an iterator
    of pileup file lines. Each line
    parsed and returned as a Locus 
    object."""

    def __init__(self, filename=None):
        if filename:
            self.infile = open(filename, 'r')
        else:
            self.infile = sys.stdin

    def __iter__(self):
        self.infile.seek(0)
        return self

    def __next__(self):
        # returns the next locus
        line = self.infile.readline()
        line = line.replace('\n','')
        line = line.split('\t')
        chrom   = line[0]
        coord   = line[1]
        ref_base = line[2]
        sample_specific_data = line[3:]
        return Locus(chrom, coord, ref_base, sample_specific_data)

class Amplification_matrix:
    """Either generates or loads from file
    a matrix of probabilities that any allele 
    will be amplified from any genotype. This
    excludes the effect of allelic dropout"""

    def __init__(self, filename=None, fp_error=0):
        # First two indeces are genotype (symmetric), third is intermediate allele
        self.matrix = np.zeros(4,4,4)
        if filename == None:
            for i in range(4):
                for j in range(4):
                    for k in range(4):
                        if i == j :
                            self.matrix[i,j,k] = 1 - fp_error if i == k else fp_error/3
                        else:
                            self.matrix[i,j,k] = 0.5 - fp_error/3 if i == k or j == k else fp_error/3
        
        else:
            f = open(filename, 'r')
            for k in range(4):
                line = f.readline().split(',')
                for col in range(len(line)):
                    i = np.floor(col/4)
                    j = col % 4
                    self.matrix[i,j,k] = line[col]
            f.close()

class VCF_file:
    """ Reads VCF file. Can be used to read
    germline mutations. """

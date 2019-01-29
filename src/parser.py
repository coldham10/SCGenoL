import sys
from locus import Locus
from collections import Iterator
import numpy as np

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
        #TODO: possibly unpredictable behaviour with stdin
        self.infile.seek(0)
        return self

    def __next__(self):
        # returns the next locus
        try:
            line = self.infile.readline()
            line = line.replace('\n','')
            line = line.split('\t')
            chrom    = line[0]
            coord    = line[1]
            ref_base = line[2]
            read_depths = np.array(line[3::3])
            reads = np.array([list(cell_reads) for cell_reads in line[4::3]])
            quals = np.array([list(cell_quals) for cell_quals in line[5::3]])

        except IndexError:
            raise StopIteration
        return Locus(chrom, coord, ref_base, read_depths, reads, quals)

class Amplification_matrix:
    """Either generates or loads from file
    a matrix of probabilities that any allele 
    will be amplified from any genotype. This
    excludes the effect of allelic dropout"""

    def __init__(self, filename=None, fp_error=0):
        # First two indices are genotype (symmetric), third is intermediate allele
        self.matrix = np.zeros((4,4,4))
        if filename == None:
            for i in range(4):
                for j in range(4):
                    for k in range(4):
                        if i == j :
                            self.matrix[i,j,k] = 1 - fp_error if i == k else fp_error/3
                        else:
                            self.matrix[i,j,k] = 0.5 - fp_error/3 if i == k or j == k else fp_error/3
        
        else:
        #TODO: test this reading ability
            f = open(filename, 'r')
            for k in range(4):
                line = f.readline().split(',')
                for col in range(len(line)):
                    i = np.floor(col/4)
                    j = col % 4
                    self.matrix[i,j,k] = line[col]
            f.close()

class VCF_file:
    #TODO test this unit
    """ Reads VCF file. Can be used to read
    germline mutations. """

    def __init__(self, filename=None):
        if filename:
            self.infile = open(filename, 'r')
        else:
            raise Exception("No VCF filename provided")

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

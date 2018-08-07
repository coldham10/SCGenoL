import sys
from locus import Locus
from collections import Iterator

class Pileup_file(Iterator):
    """Class implements an iterator
    of pileup file lines. Each line
    parsed and returned as a Locus 
    object."""

    def __init__(self, filename = None):
        if filename:
            self.infile = open(filename, 'r')
        else:
            self.infile = sys.stdin

    def __iter__(self):
        self.infile.seek(0)
        return self

    def next(self):
        # returns the next locus
        line = self.infile.readline()
        line = line.replace('\n','')
        line = line.split('\t')
        chrom   = line[0]
        coord   = line[1]
        ref_base = line[2]
        sample_specific_data = line[3:]
        return Locus(chrom, coord, ref_base, sample_specific_data)

#!/usr/bin/env python3

import sys
from collections import Iterator

class Cell_Locus:
    """Class holds information about a 
    particular cell at a particular locus"""

    def __init__(self, ref_base, n_reads, reads_string, quals_string):
        self.ref_base = ref_base
        self.read_depth = n_reads
        self.reads = reads_string
        self.quals = quals_string

    def link(self, next_cell):
        self.next_cell = next_cell


class Locus:
    """Class holds all information for a
    given locus on the reference, from 
    the pileup file"""

    def __init__(self, chrom, coord, ref_base, pileup_data):
        self.chrom = chrom
        self.coord = coord
        self.ref_base = ref_base
        self.parse_cell_data(pileup_data)
        
    def parse_cell_data(self, data):
        # Pileup data has three fields per cell
        data_bycell = [(data[3*i], data[3*i+1], data[3*i+2]) for i in range(int(len(data)/3))]
        first = data_bycell.pop(0)
        self.first_cell = Cell_Locus(self.ref_base, first[0], first[1], first[2])
        current_cell = self.first_cell
        # Create a linked list of cell data at this locus
        for cell in data_bycell:
            next_cell = Cell_Locus(self.ref_base, cell[0], cell[1], cell[2])
            current_cell.link(next_cell)
            current_cell = next_cell
            
        

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

#!/usr/bin/env python3

import sys
from collections import Iterator
import numpy as np

class Cell:
    """Class holds information about a 
    particular cell at a particular locus"""

    def __init__(self, ref_base, n_reads, reads_string, quals_string):
        self.ref_base = ref_base
        self.read_depth = n_reads
        self.reads = reads_string
        self.quals = quals_string
        #TODO: read cleaning like monovar, e.g. remove start and end symbols, change dots to ref
        #TODO: convert reads to numeric: ACGT -> 0123

    def link(self, next_cell):
        self.next_cell = next_cell

    def l_genotype_from_reads(self, g, amp_p_mat, seq_p_mat, p_ado):
        #TODO: test this method
        #NB: this is the monovar algorithm if seq_p_mat[a][b] = e/3 for a != b
        """Calculates the log likelihood of the read 
        data conditional on the cell genotype being g
        (g = 0,1 or 2 variant alleles) """
        if g == 1:
            return self.l_genotype_het(amp_p_mat, seq_p_mat, p_ado)
        bases = [0,1,2,3]
        ref = self.ref_base
        alt_bases = [0,1,2,3].remove(ref)
        # Likelihoods by underlying genotype.
        # Note that only g=0 specifies a unique genotype
        log_likelihoods = np.zeros(4)
        for i in range(self.read_depth):
            read = self.reads[i]
            if g == 0:
                # Summing over intermediate alleles to account for amplification errors
                ls_read = [amp_p_mat[ref,ref,b] * seq_p_mat[b,read] for b in bases]
                log_likelihoods[ref] += np.log(sum(ls_read))
            elif g == 2:
                # Likelihoods calculated for all possible homozygous variant genotypes
                for alt in alt_bases:
                    # Summing over intermediate alleles to account for amplification errors
                    ls_read = [amp_p_mat[alt, alt, b] * seq_p_mat[b, read] for b in bases]
                    log_likelihoods[alt] += np.log(sum(ls_read))
        return np.logaddexp.reduce(log_likelihoods)


    def l_genotype_het(self, amp_p_mat, seq_p_mat, p_ado):
        """ Calculates the log likelihood of the reads
        at this locus for this cell conditional on the cell
        being heterozygous reference/variant. Note that by
        the infinite sites model the heterozygous variant/variant
        case is ignored """
        
        bases = [0,1,2,3]
        ref = self.ref_base
        alt_bases = [0,1,2,3].remove(ref)
        # We keep track of the likelihoods for each possible genotype given allelic dropout or not
        l_ref_dropped = np.zeros(4)
        l_alt_dropped = 0
        l_no_ado = np.zeros(4)
        for i in range(self.read_depth):
            read = self.reads[i]
            ls_homo_ref = [amp_p_mat[ref, ref, b] * seq_p_mat[b, read] for b in bases]
            l_alt_dropped += np.log(sum(ls_alt_dropped))
            for alt in alt_bases:
                # amp_p_mat is the probability that the intermediate allele is amplified and fed to the sequencer
                # this already accounts for the reduced likelihood of either allele being chosen compared to the homozygous case
                ls_hetero = [amp_p_mat[alt, ref, b] * seq_p_mat[b, read] for b in bases]
                l_no_ado[alt] += np.log(sum(ls_hetero))

                ls_homo = [amp_p_mat[alt, alt, b] * seq_p_mat[b, read] for b in bases]
                l_ref_dropped[alt] += np.log(sum(ls_homo))

        # An allelic dropout of ref or alt are equally likely:
        l_ref_dropped   = np.add(l_ref_dropped, np.log(0.5))
        l_alt_dropped  += np.log(0.5)
        l_ado           = np.logaddexp(l_ref_dropped, l_alt_dropped)
        l_ado           = np.add(l_ado, np.log(p_ado))
        l_no_ado        = np.add(l_no_ado, np.log(1-p_ado))
        l_by_alt        = np.logaddexp(l_ado, l_no_ado)
        return np.logaddexp.reduce(l_by_alt)


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
        self.first_cell = Cell(self.ref_base, first[0], first[1], first[2])
        current_cell = self.first_cell
        #TODO: is this the best way?
        # Create a linked list of cell data at this locus
        for cell in data_bycell:
            next_cell = Cell(self.ref_base, cell[0], cell[1], cell[2])
            current_cell.link(next_cell)
            current_cell = next_cell

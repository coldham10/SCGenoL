import sys
from collections import Iterator
import numpy as np
import utils

class Cell:
    """Class holds information about a 
   particular cell at a particular locus"""

    def __init__(self, ref_base, n_reads, reads_string, quals_string):
        self.ref_base = ref_base
        self.read_depth = int(n_reads)
        self.reads = reads_string
        self.quals = quals_string
        utils.convert_pileup_notation(self)

    def l_genotype_from_SC_reads(self, g, amp_p_mat, p_ado):
        """Calculates the log likelihood of the read 
        data conditional on the cell genotype being g
        (g = 0,1 or 2 variant alleles) """
        if g == 1:
            return self.l_genotype_het(amp_p_mat, p_ado)
        bases = [0,1,2,3]
        ref = self.ref_base
        alt_bases = [0,1,2,3]
        alt_bases.remove(ref)
        # Likelihoods by underlying genotype.
        # Note that only g=0 specifies a unique genotype
        log_likelihoods = np.zeros(4)
        for i in range(self.read_depth):
            read = self.reads[i]
            qual = self.quals[i]
            # Probability that intermediate allele g is sequenced as h
            seq_p = lambda g, h : qual if g == h else (1 - qual) /3
            if g == 0:
                # Summing over intermediate alleles to account for amplification errors
                ls_read = [amp_p_mat[ref,ref,b] * seq_p(b,read) for b in bases]
                call_likelihoods = np.zeros(4)
                call_likelihoods[ref] = sum(ls_read)
                log_likelihoods += np.log(call_likelihoods)
            elif g == 2:
                call_likelihoods = np.zeros(4)
                # Likelihoods calculated for all possible homozygous variant genotypes
                for alt in alt_bases:
                    # Summing over intermediate alleles to account for amplification errors
                    ls_read = [amp_p_mat[alt, alt, b] * seq_p(b, read) for b in bases]
                    call_likelihoods[alt] = sum(ls_read)
                log_likelihoods += np.log(call_likelihoods)
        
        return np.logaddexp.reduce(log_likelihoods)


    def l_genotype_het(self, amp_p_mat, p_ado):
    #TODO: test this unit
        """ Calculates the log likelihood of the reads
        at this locus for this cell conditional on the cell
        being heterozygous reference/variant. Note that by
        the infinite sites model the heterozygous variant/variant
        case is ignored """
        
        bases = [0,1,2,3]
        ref = self.ref_base
        alt_bases = [0,1,2,3]
        alt_bases.remove(ref)
        # We keep track of the likelihoods for each possible genotype given allelic dropout or not
        l_ref_dropped = np.zeros(4)
        l_alt_dropped = 0
        l_no_ado = np.zeros(4)
        for i in range(self.read_depth):
            read = self.reads[i]
            qual = self.quals[i]
            seq_p = lambda g, h : qual if g == h else (1 - qual) /3
            # Case: alt allele dropped
            ls_homo_ref = [amp_p_mat[ref, ref, b] * seq_p(b, read) for b in bases]
            l_alt_dropped += np.log(sum(ls_homo_ref))
            l_call_ado = np.zeros(4)
            l_call_noado = np.zeros(4)
            for alt in alt_bases:
                # Case: ref allele dropped
                ls_homo = [amp_p_mat[alt, alt, b] * seq_p(b, read) for b in bases]
                l_call_ado[alt] = sum(ls_homo)
                # Case: no ado event
                ls_hetero = [amp_p_mat[alt, ref, b] * seq_p(b, read) for b in bases]
                l_call_noado[alt] = sum(ls_hetero)
            l_ref_dropped += np.log(l_call_ado)
            l_no_ado += np.log(l_call_noado)

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

    def __init__(self, chrom, coord, ref_base, pileup_data, germline_data=None):
        self.chrom = chrom
        self.coord = coord
        self.ref_base = ref_base
        self.parse_germline_data(germline_data)
        self.parse_cell_data(pileup_data)
        self.cells = []
        
    def parse_cell_data(self, data):
        # Pileup data has three fields per cell
        data_bycell = [(data[3*i], data[3*i+1], data[3*i+2]) for i in range(int(len(data)/3))]
        for cell_data in data_bycell:
            cell = Cell(self.germ_ref_base, cell_data[0], cell[1], cell[2])
            cells.append(cell)


    def parse_germline_data(self, data):
        #TODO: testing needed
        if data == None:
            self.germ_ref_base = self.ref_base
            self.germ_SNV = None
        else:
            (alt, g) = data
            self.germ_ref_base = alt if g == 2 else self.ref_base
            self.germ_SNV = alt


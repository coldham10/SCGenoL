#!/usr/bin/env python3

from parser import Pileup_file

#recalibrate qualities?  nielsen & SCCaller
#SFS ('individuals' could be cells)? see nielsen et al 2012 (SNP Calling, Genotype Calling, and Sample Allele...)

# monovar method


if __name__ == "__main__":
    p = Pileup_file()
    next(p)
    l = next(p)
    print(l.first_cell.next_cell.quals)

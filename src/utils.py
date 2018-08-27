""" Various functions for the genL program, including
read conversion..."""

import re
import numpy as np

base_dict =  {'a': 0,
              'A': 0,
              'c': 1,
              'C': 1,
              'g': 2,
              'G': 2,
              't': 3,
              'T': 3,
              'n': -1,
              'N': -1}

def convert_pileup_notation(cell):
    """ Cleans a cell locus of read start and end markers and
    converts all read information into numerical form"""
    cell.ref_base = base_dict[cell.ref_base]
    if cell.read_depth == 0:
        cell.reads = []
        cell.quals = []
        return
    clear_indels(cell)
    reads_to_numeric(cell)
    decode_quals(cell)
    clear_ambigs(cell)

def clear_ambigs(cell):
    """Removes ambiguous 'N' reads and their quals"""
    n_pos = [i for i in range(len(cell.reads)) if cell.reads[i] < 0]
    for pos in n_pos:
        del cell.reads[pos]
        del cell.quals[pos]
    cell.read_depth = len(cell.reads)

def clear_indels(cell):
    #TODO: handle indel quals. (?)
    """Removes indels from the read string of a cell"""
    ins_re = r"(\+[0-9]+)"
    del_re = r"(-[0-9]+)"
    insert_list = re.split(ins_re, cell.reads)
    cleaned = insert_list[0::2]
    sizes = insert_list[1::2]
    for i in range(len(sizes)):
        size = int(sizes[i][1:])
        cleaned[i+1] = cleaned[i+1][size:]
    cleaned = "".join(cleaned)
    delete_list = re.split(del_re, cleaned)
    cleaned = delete_list[0::2]
    sizes = delete_list[1::2]
    for i in range(len(sizes)):
        size = int(sizes[i][1:])
        cleaned[i+1] = cleaned[i+1][size:]
    cleaned = "".join(cleaned)
    cell.reads = cleaned

def reads_to_numeric(cell):
    """Converts the pileup read string to numeric.
    ACGT -> 0123 respectively"""
    pileup_dict =  base_dict.copy()
    pileup_dict['.'] = pileup_dict[','] = pileup_dict['*'] = cell.ref_base
    cell.reads = [ pileup_dict[r] for r in cell.reads if r in pileup_dict ]

def decode_quals(cell, qual_offset=33):
    phred = [ord(q)-qual_offset for q in cell.quals]
    p_correct = [1-np.exp(-q/10) for q in phred]
    cell.quals = p_correct

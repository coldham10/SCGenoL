""" Various functions for the genL program, including
read conversion..."""

import re

def clear_indels(cell):
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

def convert_pileup_notation(cell):
    """ Cleans a cell locus of read start and end markers and
    converts all read information into numerical form"""
    pileup_dict =  {'a': 0,
                    'A': 0,
                    'c': 1,
                    'C': 1,
                    'g': 2,
                    'G': 2,
                    't': 3,
                    'T': 3}
    cell.ref_base = pileup_dict[cell.ref_base]
    pileup_dict['.'] = pileup_dict[','] = pileup_dict['*'] = cell.ref_base


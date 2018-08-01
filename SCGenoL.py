#!/usr/bin/env python3

from pileup import Pileup_file


if __name__ == "__main__":
    p = Pileup_file()
    next(p)
    l = next(p)
    print(l.first_cell.next_cell.quals)

#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import sys
import pandas as pd
from Bio import AlignIO

def main():
    file = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    alignment = AlignIO.read(open(file), 'fasta')
    for record in alignment :
        print('>'+record.id)
        print(record.seq[start:end])

if __name__ == '__main__':
    main()

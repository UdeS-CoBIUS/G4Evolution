# -*- coding: utf-8 -*-:v

import os
import sys
import argparse
import subprocess
import multiprocessing
from functools import partial

path = sys.argv[1]

cmd = "ls FastaOrtho"
listOrtho = subprocess.check_output(cmd, shell=True)
listOrtho = listOrtho.decode("utf-8")
listOrtho = listOrtho.split('\n')
if listOrtho[-1] == '':
    listOrtho = listOrtho[:-1]

def kAlign(oN):
    cmd = "/home/vana2406/software/kalign/bin/kalign -i ~/scratch/pG4Evolution/FastaOrtho/"+oN+" -o ~/scratch/pG4Evolution/AlignmentOrtho/"+oN+" -f fasta"
    subprocess.check_output(cmd, shell=True)
    print(oN)

p = multiprocessing.Pool(31)
func = partial(kAlign)
data = p.map(func, [ i for i in listOrtho ])
p.close()

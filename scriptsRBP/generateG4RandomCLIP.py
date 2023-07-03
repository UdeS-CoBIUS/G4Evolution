#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pprint import pprint

spList = ['Hsap', 'Mmus', 'Cele', 'Drer', 'Dmel', 'Scer']

dicoStat = {'Sp': [], 'idTr': []}
for sp in spList:
    print(sp)
    pG4df = pd.read_csv('/home/anais/Documents/Projet/G4RBPevo/data/sp/'+sp+'_GAIA_Random.csv', sep=',', names=['trID', 'trStart', 'trEnd', 'chromosome', 'strand', 'pG4Start', 'pG4End', 'Sp'], index_col=None)
    pG4df['pG4Start'] = pG4df['pG4Start'].fillna(0)
    pG4df['pG4End'] = pG4df['pG4End'].fillna(0)
    pG4df['LengthTr'] = pG4df['trEnd'] - pG4df['trStart'] +1
    pG4df['LengthpG4'] = pG4df['pG4End'] - pG4df['pG4Start'] +1

    groups = pG4df.groupby('trID')
    for name, group in groups:
        dicoStat['Sp'].append(sp)
        dicoStat['idTr'].append(name)

dfStat = pd.DataFrame(data=dicoStat)

#first number = nb pG4 in exons and second number = pG4 in introns
G4Sp = {'Cele': 610+1858, 'Dmel': 3317+29758,
        'Drer': 4455+57477, 'Hsap': 54679+1253928,
        'Mmus': 29261+672515, 'Scer': 26}

dicoSpnbG4Tr = {}
for sp in spList:
    nbTr = len(list(set( dfStat[dfStat.Sp == sp]['idTr'] )))
    listnbG4Tr = [0]*nbTr
    nbG4 = G4Sp[sp]
    while nbG4 != 0:
        rd = random.randint(0, nbTr-1)
        listnbG4Tr[rd] += 1
        nbG4 -= 1
    dicoSpnbG4Tr[sp] = listnbG4Tr

for randomNb in [1,2,3,4,5,6,7,8,9,10]:
    print(randomNb)
    for sp in spList:
        listlineOutputFile = []
        dicoSpRandom = {'Start': [], 'End': [], 'Chromosome': [], 'Strand': []}
        print(sp)
        pG4df = pd.read_csv('/home/anais/Documents/Projet/G4RBPevo/data/sp/'+sp+'_GAIA_Random.csv', sep=',', names=['trID', 'trStart', 'trEnd', 'chromosome', 'strand', 'pG4Start', 'pG4End', 'Sp'], index_col=None)
        del pG4df['pG4Start']
        del pG4df['pG4End']
        pG4df = pG4df.drop_duplicates(keep='last')
        pG4df = pG4df.reset_index()

        if sp in ['Scer', 'Cele', 'Dmel']:
            lengthG4 = 20
        elif sp == 'Drer':
            lengthG4 = 30
        else:
            lengthG4 = 60

        for index, row in pG4df.iterrows():
            for i in range(0,dicoSpnbG4Tr[sp][index]):
                rd = random.randint(row.trStart, row.trEnd-5)
                try:
                    listlineOutputFile.append('chr'+str(row.chromosome)+'\t'+str(rd)+'\t'+str(rd+lengthG4)+'\t.\t.\t'+row.strand)
                except:
                    print(str(row.chromosome))
                    print(str(rd))
                    print(str(rd+lengthG4))
                    print(str(row.strand))
                    print()
                    print('chr'+str(row.chromosome)+'\t'+str(rd)+'\t'+str(rd+lengthG4)+'\t.\t.\t'+row.strand)


        output = open('/home/anais/Documents/Projet/G4RBPevo/data/Random/'+sp+'_pG4coord_'+str(randomNb)+'.bed', "a")
        output.write('\n'.join(listlineOutputFile))
        output.close()

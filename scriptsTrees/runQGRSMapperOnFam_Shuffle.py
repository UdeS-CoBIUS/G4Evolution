#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import sys
import subprocess
import pandas as pd
import multiprocessing
from Bio import AlignIO
from pprint import pprint
from functools import partial

def main(path, t):

    t = t.split('.')[0]

    res = {'idRegion': [], 'G4NN': [], 'G4H': [], 'cGcC': [], 'qgrs': [], 'gene': []}
    try:
        alignment = AlignIO.read(open(path+'AlignmentOrthoShuffle/'+t+'.fa'), 'fasta')
        # alignment = AlignIO.read(open(path+'AlignmentOrtho/Euk/Alignment_11076.fa'), 'fasta')
    except:
        print('No alignment for this tree', path+'AlignmentOrthoShuffle/'+t+'.fa')
    else:
        with open(path+'mergedpG4FamiliesPecentrOrtho/'+t+'.csv') as f:
            content = f.read()
            for l in content.split('\n'):
                l = l.split('\t')
                if l[0] != 'Family' and l[0] != '':
                    start = int(l[1])
                    end = int(l[2])

                    for record in alignment:
                        res['idRegion'].append(t+':'+str(start)+'-'+str(end))
                        res['gene'].append(record.id)
                        res['qgrs'].append(0)
                        res['G4NN'].append(0)
                        res['G4H'].append(0)
                        res['cGcC'].append(0)

                        if str(record.seq[start:end]).translate(None, '-') != '' and len(str(record.seq[start:end]).translate(None, '-')) >= 20:
                            if t+':'+str(start)+'-'+str(end) == '424:2679-2946':
                                print('>'+record.id)
                                print(str(record.seq[start-1:end]).translate(None, '-'))
                            output = open(path+'tmp.fa', "w")
                            output.write('>'+record.id+'\n')
                            output.write(str(record.seq[start-1:end]).translate(None, '-')+'\n')
                            output.close()

                            qrgs = subprocess.check_output("/home/vana2406/software/qgrs-cpp/qgrs -i "+path+"tmp.fa", shell=True)
                            # qrgs = subprocess.check_output("/home/vana2406/software/qgrs-cpp/qgrs -i "+path+"tmp.fa", shell=True)
                            qrgs = qrgs.decode("utf-8")
                            qrgsRes = qrgs.split('\n')
                            try:
                                qrgsRes.remove('ID          T1  T2  T3  T4 TS  GS  SEQ')
                            except:
                                try:
                                    qrgsRes.remove('ID         T1 T2 T3 T4 TS  GS  SEQ')
                                except:
                                    try:
                                        qrgsRes.remove('ID        T1T2T3T4 TS  GS  SEQ')
                                    except:
                                        try:
                                            qrgsRes.remove('ID           T1   T2   T3   T4 TS  GS  SEQ')
                                        except:
                                            print(qrgsRes)

                            qrgsRes.remove('')
                            qrgsRes.remove('--------------------------------------------------------------------------------------------')

                            if len(qrgsRes)-1 > 0:
                                res['qgrs'][-1] = 1

                            g4rnascreener = subprocess.check_output("/home/vana2406/software/"+\
                            "g4rna_screener/screen.py "+path+"tmp.fa -a /home/vana2406/"+\
                            "software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 "+\
                            "-c description cGcC G4H G4NN sequence start end -e", shell=True)
                            g4rnascreener = g4rnascreener.decode("utf-8")
                            g4rnascreener = g4rnascreener.split('\n')
                            for w in g4rnascreener:
                                w = w.split('\t')
                                if w[0] != '':
                                    if 'FBgn0003345' in w[1] and start == 36725 and end == 36988:
                                        print(w)
                                    # print(w)
                                    if float(w[2]) >= 4.5:
                                        res['cGcC'][-1] = 1
                                    if float(w[3]) >= 0.9:
                                        res['G4H'][-1] = 1
                                    if float(w[7]) >= 0.5:
                                        res['G4NN'][-1] = 1
        tmp = pd.DataFrame(data=res).drop_duplicates(subset=None, keep='first', inplace=False)
        tmp.to_csv(path_or_buf=path+'FamG4PredShuffle/'+t+'.csv', header=True, index=None, sep='\t')

if __name__ == '__main__':
    path = sys.argv[1]

    cmd = "ls FastaOrthoShuffle"
    listOrtho = subprocess.check_output(cmd, shell=True)
    listOrtho = listOrtho.decode("utf-8")
    listOrtho = listOrtho.split('\n')
    if listOrtho[-1] == '':
        listOrtho = listOrtho[:-1]

    p = multiprocessing.Pool(15)
    func = partial(main, path)
    data = p.map(func, [ i for i in listOrtho ])
    p.close

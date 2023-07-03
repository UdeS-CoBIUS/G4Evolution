# -*- coding: utf-8 -*-:v

import os
import math
import argparse
from Bio import SeqIO

"""
Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement of computation.

Date:
	April 2022

Description:
	This script aims to retrieve and create fasta file for all orthologs groups
    created with orthologyGraph.py.

"""

def writeFasta(outDir, dicoFasta):
    """From a dictionary {id : seq}, write a fasta file.

    :param outDir: name of the directory where the output file need to be writen.
    :type outDir: string
    :param dicoFasta: {id : seq}
    :type dicoFasta: dictionary
    """
    output = open(outDir, "w")
    for id in dicoFasta:
        output.write('>' + id + "\n")
        nbLine = math.ceil( float( len(dicoFasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            seq = str(dicoFasta[id])
            output.write(str(seq[cpt1:cpt2]) + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

def getFasta(path, listGene, dicoGeneList, outDir):
    """Retrieve fasta of each gene from an ortholog group.

    :param path: name of the main directory.
    :type path: string
    :param listGene: list of genes from a group of orthologs.
    :type listGene: list
    :param dicoGeneList: {gene : sp}.
    :type dicoGeneList: dictionary
    :param outDir: name of the output fasta file.
    :type outDir: string
    """
    dicoFastaTmp = {}
    for gene in listGene:
        if gene in dicoGeneList:
            sp = dicoGeneList[gene]
            filename = path + '/data/' + sp + '/Sequences_Gene_WT.fa'
            geneUnspliced = {record.id: record for record in SeqIO.parse(filename, "fasta")}
            parsedGeneUnspliced = {}
            infoGene = {}
            for seqId in geneUnspliced:
                parsedGeneUnspliced[ seqId.split(':')[0] ] = geneUnspliced[seqId]
                infoGene[ seqId.split(':')[0] ] = seqId
            if gene in parsedGeneUnspliced:
                dicoFastaTmp[ parsedGeneUnspliced[gene].description ] = parsedGeneUnspliced[gene].seq
    writeFasta(outDir, dicoFastaTmp)

def main(dicoFinal, dicoGeneListbySp, path, idFam):
    print(idFam)
    outputDir = path + 'FastaOrtho/'+idFam+'.fa'
    getFasta(path, dicoFinal[idFam], dicoGeneListbySp, outputDir)

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'getFastaTree')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-oN', '--orthoName', default = 13111)
    parser.add_argument ('-grp', '--spGrp', default = 'Arc')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    oN = arg.orthoName
    grp = arg.spGrp
    with open(path+'listOrtho/'+grp+'/'+str(oN)+'.txt') as f:
        content = f.read()
        content = content.split('\n')
    main(content, grp+'/grpOrtho_'+str(oN))

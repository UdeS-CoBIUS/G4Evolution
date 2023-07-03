#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import re
import math
import argparse
import subprocess
import pandas as pd
import multiprocessing
from Bio import AlignIO
from pprint import pprint
from functools import partial

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	May 2022

Description:
	This script aims to read one MSA and to find pG4 families. To do that,
	pG4 of genes from the MSA are retrevied. Then, we check if those pG4 are
	overlapping in the MSA. If they do, it makes a pG4 family.

"""

def readGTF(file, gene):
	"""Reads a gtf file to import information on genes in the alignment.

	:param file: name of the pG4 file.
	:type file: str
	:param gene: gene id.
	:type gene: str
	"""
	with open(file) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if gene in l and l.split('\t')[2] == 'gene':
				l = l.split('\t')
				biotype = ''
				for attribute in l[8].split(';'):
					if re.search("gene_biotype", attribute):
						biotype = attribute.split('"')[1]
				attributes = l[8].split(';')
				geneID = attributes[0].split('"')[1]
				start = min(int(l[3]), int(l[4]))
				end = max(int(l[3]), int(l[4]))
				d = {'Start' : start,
					'End' : end,
					'Strand': l[6],
					'Chr': l[0],
					'Biotype': biotype,
					'Gene' : geneID}
				return d

def readpG4(file, gtf):
	"""Improt all pG4 from a specic gene.

	:param file: name of the pG4 file.
	:type file: str
	:param gtf: gets information of genes we want to check if they have pG4.
	:type gtf: dictionary

	:returns: pG4s, dicoCoordpG4, {'CoordAlignment','Gene','CoordReal','Species','Chromosome','Strand'}
	:rtype: list, dictionary
	"""
	pG4s = []
	dicoCoordpG4 = {}
	# sp = file.split('/')[7]
	sp = file.split('/')[5]
	with open(file) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			w = l.split('\t')
			if w[0] != '':
				if w[4] == gtf['Gene']:
					d = {'Gene' : w[4],
						'pG4Seq': w[10]}
					# change coords from chromosomal to 'relative' ones
					pG4Start = min(int(w[8]), int(w[9]))
					pG4End = max(int(w[8]), int(w[9]))
					if gtf['Strand'] == '+':
						d['Start'] = pG4Start - gtf['Start'] + 1
						d['End'] = pG4End - gtf['Start'] + 1
					else:
						d['Start'] = gtf['End'] - pG4End + 1
						d['End'] = gtf['End'] - pG4Start + 1
					pG4s.append(d)
					pG4coord = str(d['Start']) +'-'+ str(d['End'])+'-'+w[4]
					dicoCoordpG4[pG4coord] = {'CoordAlignment': [0, 0],
						'Gene': w[4],
						'CoordReal': [pG4Start, pG4End],
						'Species': sp, 'Chromosome': gtf['Chr'], 'Strand': gtf['Strand']}
				# else:
				# 	print(w)
	# pprint(dicoCoordpG4)
	return pG4s, dicoCoordpG4

def shannon_entropy(list_input):
	"""
		Taken from https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py
	"""
	# Get only the unique bases in a column
	unique_base = set(list_input)
	# if '-' in unique_base:
	# 	unique_base.remove("-")
	M = len(list_input)
	entropy_list = []
	# Number of residues in column
	for base in unique_base:
		n_i = list_input.count(base) # Number of residues of type i
		P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
		entropy_i = P_i*(math.log(P_i,2))
		entropy_list.append(entropy_i)
	sh_entropy = -(sum(entropy_list))
	return sh_entropy

def computeID(list_input):
	max = 0
	nbSeq = len(list_input)
	for nt in set(list_input):
		ntPerc = float(list_input.count(nt)) / float(nbSeq)
		if ntPerc > max:
			max = ntPerc
	return max

def computeRationSeq(list_input):
	nbSeq = len(list_input)
	ntPerc = float(list_input.count('-')) / float(nbSeq)
	return 1-ntPerc

def readAlignment(alignmentFile, dicoGeneList, path, tree):
	"""Reads the alignment file and create a table.

	A table is created where each cells/element correspond to a pG4 (1),
	a non pG4 (0) or a gap (-). In order to do this, we first find genes
	information (their pG4, which gene is at which postition) and then we
	browse all alignement position and check if there is a pG4 or no. The
	parsed alignment is printed out of a dataframe into an output file.

	:param alignmentFile: name of the alignment file.
	:type alignmentFile: str
	:param dicoGeneList: {geneID : specie name}
	:type padicoGeneListth: dictionary
	:param path: main forlder path.
	:type path: str
	:param tree: id of the tree we are analyzing the alignment
	:type tree: str

	:returns: dicoCoordpG4Tot {'CoordAlignment','Gene','CoordReal','Species','Chromosome','Strand'}
	:rtype: dictionary
	"""
	pG4s = {}
	orderGene = [] # remember the order of genes in the alignment
	gtf = {}
	cpt_Alg = {} # compteur alignment
	cpt = 0
	entropy = []
	identity = []
	ratioSeq = []
	dicoCoordpG4Tot = {}
	alignment = AlignIO.read(open(alignmentFile), 'fasta')
	print('Alignement length : ', alignment.get_alignment_length())
	for record in alignment :
		gene = record.id.split(':')[0]
		#browse all genes in the alignment
		cpt_Alg[gene] = 0
		orderGene.append(gene)
		sp = dicoGeneList[gene]
		gtfFile = path + 'data/'+sp+'/'+sp+'.gtf'
		pG4File = path + 'data/'+sp+'/Gene_pG4WT.csv'
		gtf[gene] = readGTF(gtfFile, gene)
		cpt += 1
		if gene not in pG4s:
			pG4s[gene] = []
		pG4tmp, dicoCoordpG4 = readpG4(pG4File, gtf[gene])
		dicoCoordpG4Tot.update(dicoCoordpG4)
		if pG4tmp:
			for pG4 in pG4tmp:
				pG4s[gene].append(pG4)
	listCol = []
	# pprint(pG4s)
	# pprint(dicoCoordpG4Tot)
	for aliPos in range(0,alignment.get_alignment_length()):
		# browse all position in the alignment
		cptGene = 0 # reinitialisation of the gene compteur
		col = [0] * cpt
		entropy.append(shannon_entropy(list(alignment[:, aliPos])))
		identity.append(computeID( list(alignment[:, aliPos]) ))
		ratioSeq.append( computeRationSeq( list(alignment[:, aliPos]) ) )
		for g in range(0, len(alignment[:, aliPos]) ):
			#browse all gene the position i of the alignement
			if alignment[:, aliPos][g] != '-':
				cpt_Alg[ orderGene[cptGene] ] += 1
				for pg4 in pG4s[ orderGene[cptGene] ]:
					# browese all pG4 of a gene to find if there are at that
					# position
					if pg4 : # genes can have no pG4 so we need to test it
						# if orderGene[cptGene] == 'OE_3925R':
						# 	print(pg4['Start'], cpt_Alg[ orderGene[cptGene] ])
						if pg4['Start'] <= cpt_Alg[ orderGene[cptGene] ] \
							and cpt_Alg[ orderGene[cptGene] ] <= pg4['End']:
							pG4coord = str(pg4['Start'])+'-'+str(pg4['End'])+'-'+orderGene[cptGene]
							if dicoCoordpG4Tot[pG4coord]['CoordAlignment'][0] == 0:
								dicoCoordpG4Tot[pG4coord]['CoordAlignment'][0] = aliPos+1
								dicoCoordpG4Tot[pG4coord]['CoordAlignment'][1] = aliPos+1
							else:
								dicoCoordpG4Tot[pG4coord]['CoordAlignment'][1] = aliPos+1
							if orderGene[cptGene] == 'FBgn0038516':
								print('Patate')
							col[cptGene] = 1
						else:
							if col[cptGene] != 1:
								col[cptGene] = 0
					else:
						if col[cptGene] != 1:
							col[cptGene] = 0
			else:
				col[cptGene] = '-'
			cptGene += 1 # uptdate of the gene compteur
		listCol.append(col)
	# orderGene.append('Entropy')
	df = pd.DataFrame(listCol, columns=orderGene)
	dftmp = pd.DataFrame()
	dftmp = dftmp.append(df)
	dftmp['Entropy'] = entropy
	dftmp['Identiy'] = identity
	dftmp['RatioSeq'] = ratioSeq
	dftmp.to_csv(path_or_buf=path+'ProcessedAlignOrtho/'+tree.split('.')[0]+'.csv', header=True, index=None, sep='\t')
	return dicoCoordpG4Tot

def importGeneList(geniListDir):
	"""Imports all gene list and sp into a dictionary.

	:param geniListDir: name of the directory which contains each species as
		a file with all their genes.
	:type geniListDir: str

	:returns: dicoGeneListbySp {gene:sp}
	:rtype: dictionary
	"""
	dicoGeneListbySp = {}
	for path, dirs, files in os.walk(geniListDir):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = geniListDir + filename
			sp = filename.split('.')[0]
			with open(inputfile) as f:
				content = f.read()
				lines = content.split('\n')
				for l in lines:
					dicoGeneListbySp[l] = sp
	return dicoGeneListbySp

def main(path, tree):
	"""Main functions calling.

	This script have 3 main parts. First we create a dictionary of genes with
	species as values. Then, using this dictionary we could find if there is pG4
	in those genes and create a dataFrame with lines for each genes and columns
	for each position in the alignment. Each cells correspond to a position for
	gene in the alignment can be : 0 -> no pG4, 1 -> pG4, - -> gap.

	:param path: main forlder path.
	:type path: str
	:param tree: id of the tree we are analyzing the alignment
	:type tree: str
	"""
	print('Current aligmnent : '+tree)
	df = pd.DataFrame()
	alignmentFile = path + 'AlignmentOrtho/'+tree
	geniListDir = path + 'Homology/'
	dicoGeneList = importGeneList(geniListDir)
	dicoCoordpG4Tot = readAlignment(alignmentFile, dicoGeneList, path, tree)
	outputF = path + 'coordpG4/'+tree.split('.')[0]+ '.csv'
	output = open(outputF, "w")
	output.write("sp\tStrand\tChromosome\tGene\tpG4Start\tpG4End\talignmentStart\talignmentEnd\n")
	for pG4 in dicoCoordpG4Tot:
		output.write(dicoCoordpG4Tot[pG4]['Species']+'\t'+dicoCoordpG4Tot[pG4]['Strand']+\
		'\t'+dicoCoordpG4Tot[pG4]['Chromosome']+'\t'+dicoCoordpG4Tot[pG4]['Gene']+'\t'+str(dicoCoordpG4Tot[pG4]['CoordReal'][0])+'\t'+\
		str(dicoCoordpG4Tot[pG4]['CoordReal'][1])+'\t'+str(dicoCoordpG4Tot[pG4]['CoordAlignment'][0])+'\t'+\
		str(dicoCoordpG4Tot[pG4]['CoordAlignment'][1])+'\n')
	output.close()

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getpG4Familly')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path

	cmd = "ls FastaOrtho"
	listOrtho = subprocess.check_output(cmd, shell=True)
	listOrtho = listOrtho.decode("utf-8")
	listOrtho = listOrtho.split('\n')
	if listOrtho[-1] == '':
	    listOrtho = listOrtho[:-1]

	p = multiprocessing.Pool(31)
	func = partial(main, path)
	data = p.map(func, [ i for i in listOrtho ])
	p.close

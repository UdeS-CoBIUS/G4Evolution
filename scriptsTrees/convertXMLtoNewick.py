#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd
from Bio import Phylo

def main(tree, path):
	inputfile = path+"prunedTree/protein_default.tree."+tree+".phyloxml.xml"
	output = path+"newickTree/"+tree+".tre"
	Phylo.convert(inputfile, "phyloxml", output, "newick")
	cmd = "sed -i s/_//g "+output
	os.system(cmd)

	#tree = Phylo.parse(output, 'newick')
	#Phylo.draw(tree)

	inputGeneList = path+"FiltredGeneListByGrpAfterNbfilter/ArcBacEuk/gene_"+tree+".txt"
	dicoGene = {}
	with open(inputGeneList) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			l = l.replace('_', '')
			dicoGene[l] = []

	inputpG4Family = path+"mergedpG4FamiliesTestPer/"+tree+".csv"
	with open(inputpG4Family) as f:
		content = f.read()
		lines = content.split('\n')
		lines = lines[1:]
		nbFam = len(lines)

		for g in dicoGene:
			dicoGene[g] = [0]*nbFam

		for l in lines:
			w = l.split('\t')
			fam = w[0]
			if fam != '':
				intFam = int(fam.split('m')[1])
				genes = w[3].split('|')
				print(intFam)
				for g in genes:
					g = g.replace('_', '')
					dicoGene[g][intFam] += 1

		col = []
		cpt = 0
		for n in lines:
			col.append('Fam'+str(cpt))
			cpt += 1

	df = pd.DataFrame(data=dicoGene)
	df.index = col
	df.to_csv(path_or_buf=path+'/Matrices/'+tree+'.csv', header=True, index=True, sep='\t')

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'prepareBlast')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-t', '--tree', default = '390858')
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	tree = arg.tree
	main(tree, path)

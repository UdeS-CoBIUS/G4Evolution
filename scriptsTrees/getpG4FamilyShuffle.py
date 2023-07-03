#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import subprocess
import pandas as pd
import multiprocessing
from pprint import pprint
from functools import partial

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	MAy 2022

Description:
	This script aims to read processed MSA and to find pG4 families. To do that,
	pG4 of genes from the MSA are retrevied. Then, we check if those pG4 are
	overlapping in the MSA. If they do, it makes a pG4 family.
	This script is running with multiprocess.

"""

def get_index_positions(list_of_elems, element):
	"""Returns the indexes of all occurrences of give element in
	the list- listOfElements

	--Taken from StackOverflow--
	"""
	index_pos_list = []
	index_pos = 0
	while True:
		try:
			# Search for item in list from indexPos to the end of list
			index_pos = list_of_elems.index(element, index_pos)
			# Add the index position in list
			index_pos_list.append(index_pos)
			index_pos += 1
		except ValueError as e:
			break
	# print(index_pos_list)
	# print('----')
	return index_pos_list

def getFamily(df):
	"""Gets pG4 familly from the alignment dataFrame.

	All columns (position of the alignment) of the dataFrame are browsed. First,
	we check if there is a 1 in the column. If yes, we create/update the family
	and save the genes that have a pG4. Then, if there is no 1, we check if there
	are gaps (-), if there is a gap we check if the last number was a 1 which
	means we are in a pG4 family but in a gap. Finally, if we are in none of
	those cases, it means there are only 0 and no family.

	:param df: contains all information about the alignment, 1 for pG4, 0 for
		no pG4 and - for a gap.
	:type df: dataFrame

	:returns: families, pG4FamillyCoord, {idFam :{'Genes', 'Start', 'End'}}
	:rtype: dictionary, list
	"""
	familyBool = False
	families = {}
	lastData = {'nb' : [], 'cpt' : []}
	prev = {}
	lastI = 0
	pG4FamillyCoord = []
	for i in range(0, df.shape[1]):
		# browse all column of the dataFrame
		if '1' in list(df.iloc[:,i]) or 1 in list(df.iloc[:,i]):
			# if there is a pG4
			if not familyBool:
				# It's the begining of a new family
				try:
					# if len(families[idFam]['Genes'])> 1:
					pG4FamillyCoord.append([start, end])
				except:
					print('first fam, no coord yet')
				start = i
				end = i
				idFam = 'Fam'+str(len(families))
			if idFam not in families:
				# create the family
				genes = []
				lst = list(df.iloc[:,i])
				for g in get_index_positions(lst, '1'):
					genes.append( list(df.iloc[:,i].index)[g] )
				if genes == []:
					for g in get_index_positions(lst, 1):
						genes.append( list(df.iloc[:,i].index)[g] )
				families[idFam] = {'Genes' : genes, 'Start' : i+1, 'End' : i+1}
			else:
				#update the family
				end += 1
				genes = []
				lst = list(df.iloc[:,i])
				for g in get_index_positions(lst, '1'):
					genes.append( list(df.iloc[:,i].index)[g] )
				for g in genes:
					if g not in families[idFam]['Genes']:
						families[idFam]['Genes'].append(g)
				families[idFam]['End']  = i+1
			familyBool = True
			for g in get_index_positions(lst, '1'):
				prev[g] = '1'
			lastI = i
		elif '-' in list(df.iloc[:,i]):
			# try to find if the gap is in the middle of a family
			lst = list(df.iloc[:,i])
			pos = get_index_positions(lst, '-')
			# tmpPrev = prev.copy()
			prevBool = False
			for g in prev:
				if g in pos:
					if prev[g] == '1':
						prevBool = True
						# and df.iloc[g,i]
						families[idFam]['End'] = lastI+1
						break
			for g in prev:
				if g not in pos and lst[g] == '0' and prevBool == False:
					familyBool = False
					# del tmpPrev[g]
					prev = {}
					break
			# prev = tmpPrev.copy()
		else:
			# there is no family so we reinitialsed our variable
			prev = {}
			familyBool = False
	return families, pG4FamillyCoord

def writeFamilies(outputF, families):
	"""Write the output containing all famillies from an alignment.

	:param outputF: name of the output file.
	:type outputF: string
	:param families: {idFam: {'Genes', 'Start', 'End' }}
	:type families: dictionary
	"""
	output = open(outputF, "w")
	output.write("Family\tStart\tEnd\tGenes\tnbGenes\n")
	for fam in families:
		output.write(fam +"\t"+ str(families[fam]['Start']) +"\t"+ \
			str(families[fam]['End']) +"\t"+\
			'|'.join(families[fam]['Genes']) +"\t"+ \
			str(len(families[fam]['Genes'])) + "\n")
	output.close()

def getSizeBetweenFam(families):
	"""Gets the length inter-familly regions.

	:param families: {idFam: {'Genes', 'Start', 'End' }}
	:type families: dictionary

	:returns: size, contains the sizes inter-familly regions
	:rtype: list
	"""
	size = []
	for i in range(0, len(families)):
		if i != len(families)-1:
			fam = 'Fam'+str(i)
			famp = 'Fam'+str(i+1)
			size.append(families[famp]['Start']-families[fam]['End']+1)
	return size

def mergeFam(families, w):
	"""Merge famillies closer then a certain window.

	:param families: {idFam: {'Genes', 'Start', 'End' }}
	:type families: dictionary
	:param w: lenght of the window between famillies.
	:type w: int

	:returns: size, contains the sizes inter-familly regions
	:rtype: list
	"""
	dico = {}
	try:
		mergedFam = [families['Fam0']]
	except:
		print('There are no pG4 fam in this gene tree.')
	else:
		for i in range(1, len(families)):
			c = 'Fam'+str(i)
			previous = mergedFam[-1]
			current = families[c]
			if previous['End']+w > families[c]['Start']:

				previous['End'] = max(previous['End'], families[c]['End'])
				current['Start'] = previous['Start']
				current['End'] = previous['End']
				current['Genes'].extend(previous['Genes'])
				current['Genes'] = list(set(current['Genes']))

				mergedFam.append( current )

			else:
				mergedFam.append(current)

		mergedFam2 = []
		for i in range(0, len(mergedFam)):
			if i != len(mergedFam)-1:
				if mergedFam[i]['Start'] != mergedFam[i+1]['Start']:
					mergedFam2.append(mergedFam[i])
				else:
					pass
			else:
				mergedFam2.append(mergedFam[i])

		for i in range(0, len(mergedFam2)):
			c = 'Fam'+str(i)
			dico[c] = mergedFam2[i]
	return dico

def main(path, percent, w, tree):
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
	processedAlignmentFile = path + 'ProcessedAlignOrtho/'+tree
	geniListDir = path + 'Homology/'
	df = pd.read_csv(processedAlignmentFile, sep='\t')
	lengthAlignment = len(df)
	print(lengthAlignment)
	del df['Entropy']
	del df['Identiy']
	del df['RatioSeq']
	df = df.T
	families, pG4FamillyCoord = getFamily(df)
	if w == 0:
		outputF = path + 'pG4FamiliesOrtho/'+tree
	else:
		if percent:
			w = (w / 100) * lengthAlignment
			outputF = path + 'mergedpG4FamiliesPecentrOrtho/'+tree
			families = mergeFam(families, w)
		else:
			outputF = path + 'mergedpG4FamiliesOrtho/'+tree
			families = mergeFam(families, w)
	pprint(families)
	print(outputF)
	writeFamilies(outputF, families)
	size = getSizeBetweenFam(families)

	size = [str(x) for x in size]
	outputF = path+'gapFamOrtho/'+tree
	output = open(outputF, "w")
	output.write('\n'.join(size))
	output.close()

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getpG4Familly')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-w', '--window', default = 0)
	parser.add_argument ('-pc', '--percent', default = False)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	percent = arg.percent
	w = float(arg.window)

	cmd = "ls ProcessedAlignOrtho | grep 'Repro'"
	listOrtho = subprocess.check_output(cmd, shell=True)
	listOrtho = listOrtho.decode("utf-8")
	listOrtho = listOrtho.split('\n')
	if listOrtho[-1] == '':
		listOrtho = listOrtho[:-1]

	p = multiprocessing.Pool(2)
	func = partial(main, path, percent, w)
	data = p.map(func, [ i for i in listOrtho ])
	p.close

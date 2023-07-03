# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd
import getOrthoFasta
import networkx as nx
import seaborn as sns
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt

"""
Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement of computation.

Date:
	April 2022

Description:
	This script aims to create orthologs 1:1 groups from pan EnsemblCompara
	homology links.

"""

def createDicoSpGrp():
	"""Create a dictionary {sp:grp}

	:returns: dicoSp, {sp:grp}
	:rtype: dictionary
	"""
	dicoSp = {'pan_troglodytes' : 'Eukaryota',
	'homo_sapiens' : 'Eukaryota',
	'pongo_abelii' : 'Eukaryota',
	'mus_musculus' : 'Eukaryota',
	'monodelphis_domestica' : 'Eukaryota',
	'ornithorhynchus_anatinus' : 'Eukaryota',
	'anolis_carolinensis' : 'Eukaryota',
	'gallus_gallus' : 'Eukaryota',
	'danio_rerio' : 'Eukaryota',
	'gasterosteus_aculeatus' : 'Eukaryota',
	'drosophila_melanogaster' : 'Eukaryota',
	'apis_mellifera' : 'Eukaryota',
	'caenorhabditis_elegans' : 'Eukaryota',
	'neurospora_crassa' : 'Eukaryota',
	'aspergillus_nidulans' : 'Eukaryota',
	'saccharomyces_cerevisiae' : 'Eukaryota',
	'schizosaccharomyces_pombe' : 'Eukaryota',
	'dictyostelium_discoideum' : 'Eukaryota',
	'arabidopsis_thaliana' : 'Eukaryota',
	'vitis_vinifera' : 'Eukaryota',
	'solanum_lycopersicum' : 'Eukaryota',
	'oryza_sativa' : 'Eukaryota',
	'physcomitrella_patens' : 'Eukaryota',
	'chlamydomonas_reinhardtii' : 'Eukaryota',
	'leishmania_major' : 'Eukaryota',
	'methanosarcina_acetivorans_c2a' : 'Archaea',
	'halobacterium_salinarum_r1' : 'Archaea',
	'hyperthermus_butylicus_dsm_5456' : 'Archaea',
	'archaeoglobus_fulgidus_dsm_4304' : 'Archaea',
	'methanobrevibacter_smithii_atcc_35061' : 'Archaea',
	'pyrococcus_horikoshii_ot3' : 'Archaea',
	'thermoplasma_acidophilum_dsm_1728' : 'Archaea',
	'sulfolobus_solfataricus_p2' : 'Archaea',
	'pyrobaculum_aerophilum_str_im2' : 'Archaea',
	'nanoarchaeum_equitans_kin4_m' : 'Archaea',
	'candidatus_korarchaeum_cryptofilum_opf8' : 'Archaea',
	'cenarchaeum_symbiosum_a' : 'Archaea',
	'aquifex_aeolicus_vf5' : 'Bacteria',
	'mycoplasma_pneumoniae_m129' : 'Bacteria',
	'staphylococcus_aureus_subsp_aureus_n315' : 'Bacteria',
	'bacillus_subtilis_subsp_subtilis_str_168' : 'Bacteria',
	'enterococcus_faecalis_v583' : 'Bacteria',
	'streptococcus_pneumoniae_tigr4' : 'Bacteria',
	'chloroflexus_aurantiacus_j_10_fl' : 'Bacteria',
	'mycobacterium_tuberculosis_h37rv' : 'Bacteria',
	'thermus_thermophilus_hb8' : 'Bacteria',
	'chlamydia_trachomatis_d_uw_3_cx' : 'Bacteria',
	'borrelia_burgdorferi_b31' : 'Bacteria',
	'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819' : 'Bacteria',
	'myxococcus_xanthus_dk_1622' : 'Bacteria',
	'geobacter_sulfurreducens_pca' : 'Bacteria',
	'wolbachia_endosymbiont_of_drosophila_melanogaster' : 'Bacteria',
	'anaplasma_phagocytophilum_str_hz' : 'Bacteria',
	'brucella_abortus_bv_1_str_9_941' : 'Bacteria',
	'neisseria_meningitidis_z2491' : 'Bacteria',
	'legionella_pneumophila_str_paris' : 'Bacteria',
	'francisella_tularensis_subsp_tularensis_schu_s4' : 'Bacteria',
	'vibrio_cholerae_o1_biovar_el_tor_str_n16961' : 'Bacteria',
	'haemophilus_influenzae_rd_kw20' : 'Bacteria',
	'yersinia_pestis_biovar_microtus_str_91001' : 'Bacteria',
	'escherichia_coli_str_k_12_substr_mg1655' : 'Bacteria'}
	return dicoSp

def importHomology(filename, dicoSp):
	"""Import orthology information into a graphe.

	The file with all orthology link, and for each pair of gene, we check if
	they both are from a species in our study, if yes and that they are not
	already inside the graphe, we add them with, then we add the edge between the
	two gene.

	:param filename: name of the file with all orthology 1:1 links.
	:type filename: string
	:param dicoSp: {sp:grp}.
	:type dicoSp: dictionary
	:returns: G, graph of orthologs built with genes as nodes and orthology links
	as edges
	:rtype: network
	"""
	G = nx.Graph()
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if l:
				words = l.split('\t')
				gene1 = words[0]
				gene2 = words[5]
				homologyT = words[4]
				specie1 = words[2]
				specie2 = words[7]
				if specie1 in dicoSp and specie2 in dicoSp:
					#if both species are among those we are studying
					if gene1 not in G:
						G.add_node(gene1, specie=specie1)
					if gene2 not in G:
						G.add_node(gene2, specie=specie2)
					G.add_edge(gene1, gene2)
	return G

def importGeneList(geniListDir):
	"""Get a dico to have for each gene their sp.

	This function aims browse all genes of all species in order to know from
	which species are a gene.

	:param geniListDir: name of the directory containing all species file.
	:type geniListDir: string
	:returns: dicoGeneListbySp, {geneId:sp}
	:rtype: dictionary
	"""
	dicoGeneListbySp = {}
	for path, dirs, files in os.walk(geniListDir):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = geniListDir + filename
			if 'txt' in inputfile:
				sp = filename.split('.')[0]
				with open(inputfile) as f:
					content = f.read()
					lines = content.split('\n')
					for l in lines:
						dicoGeneListbySp[l] = sp
	return dicoGeneListbySp

def getOrthoGrp(G, cc):
	"""Retrieve connex componant (orthologs grp) from the graph.

	All cc are retrieved out of the main graph, then they are all browsed and only
	group of 3 or more genes are kept. This threshold will be updated latter
	after the analyse of orthologs groups size. A random identifier is attributed
	for each orthologs group.

    :param G: graph of orthologs built with genes as nodes and orthology links
	as edges.
    :type G: graph
	:param cc: connex componant of the main graph.
    :type cc: subgraph
    :returns: dicoOrthoFam, {orthoID:list of orthologs}
    :rtype: dictionary
    """
	cpt = 0
	dicoOrthoFam = {}
	for familly in cc:
		famG = nx.Graph(G.subgraph(familly))
		if len(famG) >= 3:
			dicoOrthoFam['grpOrtho_'+str(cpt)] = list(famG.nodes)
			cpt += 1
	return dicoOrthoFam

def getGrpSpOrthoGrp(geneList, dicoGrpSp, dicoSpGene):
	"""Retrieve for a group of orthologs the living kingdom each genes are from.

	This function allow us to find for a group of orthologs, all the species
	group from the genes inside that ortholog group. It will allow us to check
	the number of genes for each sub grp and be more faire on choosing thresholds.

	:param geneList: list of gene from the ortholog group.
	:type geneList: list
	:param dicoGrpSp: {sp : group}.
	:type dicoGrpSp: dictionary
	:param dicoSpGene: {gene : sp}.
	:type dicoSpGene: dictionary

	:returns: name of the species group that genes are in
    :rtype: string
	"""
	dicoTree = {}
	dicoGrp = {'Archaea': [], 'Bacteria': [], 'Eukaryota': []}

	for g in geneList:
		if g in dicoSpGene:
			sp = dicoSpGene[g]
			grp = dicoGrpSp[sp]
			dicoGrp[grp].append(sp)
	for grp in dicoGrp:
		dicoGrp[grp] = len(list(set(dicoGrp[grp])))

	if dicoGrp['Archaea'] > 0 and dicoGrp['Bacteria'] > 0 and dicoGrp['Eukaryota'] > 0:
		return 'ArcBacEuk'
	elif dicoGrp['Archaea'] > 0 and dicoGrp['Bacteria'] > 0 and dicoGrp['Eukaryota'] == 0:
		return 'ArcBac'
	elif dicoGrp['Archaea'] == 0 and dicoGrp['Bacteria'] > 0 and dicoGrp['Eukaryota'] > 0:
		return 'BacEuk'
	elif dicoGrp['Archaea'] > 0 and dicoGrp['Bacteria'] == 0 and dicoGrp['Eukaryota'] > 0:
		return 'ArcEuk'
	elif dicoGrp['Archaea'] > 0 and dicoGrp['Bacteria'] == 0 and dicoGrp['Eukaryota'] == 0:
		return 'Arc'
	elif dicoGrp['Archaea'] == 0 and dicoGrp['Bacteria'] > 0 and dicoGrp['Eukaryota'] == 0:
		return 'Bac'
	elif dicoGrp['Archaea'] == 0 and dicoGrp['Bacteria'] == 0 and dicoGrp['Eukaryota'] > 0:
		return 'Euk'

def main(filename):
	"""Main function.

	A graph with orthologs 1:1 is created, from which we can then retrieve
	orthologs groups by isolationg connex componant. This script was launch a
	first time to check the size of orthology groups and then a second time
	retrieve only groups within this size. There is a different threshold because
	all groups are not equal.

	:param filename: sequence to shuffle.
	:type filename: string
	"""
	path = '/'.join(filename.split('/')[0:-2])+'/'
	dicoGrpSp = createDicoSpGrp()
	G = importHomology(filename, dicoGrpSp)
	cc = nx.connected_components(G)
	dicoOrthoFam = getOrthoGrp(G, cc)
	dicoGeneListbySp = importGeneList(path+"Homology/")
	dicoGrpLenFam = {'Grp': [], 'Len':[]}
	dicoFinal = {}

	for oN in dicoOrthoFam:
		dicoOrthoFam[oN] = list(set(dicoOrthoFam[oN]))
		grp = getGrpSpOrthoGrp(dicoOrthoFam[oN], dicoGrpSp, dicoGeneListbySp)
		dicoGrpLenFam['Len'].append(len(dicoOrthoFam[oN]))
		dicoGrpLenFam['Grp'].append(grp)
		if grp == 'Arc' and len(dicoOrthoFam[oN]) >= 7:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'ArcBac' and len(dicoOrthoFam[oN]) >= 10:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'ArcEuk' and len(dicoOrthoFam[oN]) >= 15:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'ArcBacEuk' and len(dicoOrthoFam[oN]) >= 20:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'BacEuk' and len(dicoOrthoFam[oN]) >= 15:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'Bac' and len(dicoOrthoFam[oN]) >= 10:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]
		elif grp == 'Euk' and len(dicoOrthoFam[oN]) >= 12:
			output = open(path+"listGrpOrtho/"+oN+".txt", "a")
			output.write('\n'.join(dicoOrthoFam[oN]))
			output.close()
			dicoFinal[oN] = dicoOrthoFam[oN]


	p = multiprocessing.Pool(31)
	func = partial(getOrthoFasta.main, dicoFinal, dicoGeneListbySp, path)
	data = p.map(func, [ i for i in dicoFinal ])
	p.close()

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'prepareBlast')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	filename = path+"Homology/AllSpOrto.tsv"
	main(filename)

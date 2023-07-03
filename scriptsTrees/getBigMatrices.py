#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import subprocess
import pandas as pd
from pprint import pprint

def importGeneList(geniListDir, listSp):
	"""Imports all gene list and sp into a dictionary.

	:param geniListDir: name of the directory which contains each species as
		a file with all their genes.
	:type geniListDir: str
	"""
	dicoGeneListbySp = {}
	for path, dirs, files in os.walk(geniListDir):
		# for each element of the directory to passed
		for filename in files: # for each files
			if 'txt' in filename:
				inputfile = geniListDir + filename
				sp = filename.split('.')[0]
				if sp in listSp:
					with open(inputfile) as f:
						content = f.read()
						lines = content.split('\n')
						for l in lines:
							dicoGeneListbySp[l] = sp
	return dicoGeneListbySp

def getSpParse():
	list = ['bacillus_subtilis_subsp_subtilis_str_168', 'francisella_tularensis_subsp_tularensis_schu_s4', 'leishmania_major', 'pongo_abelii', 'saccharomyces_cerevisiae', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'anaplasma_phagocytophilum_str_hz', 'pyrococcus_horikoshii_ot3', 'chlamydomonas_reinhardtii', 'vitis_vinifera', 'mus_musculus', 'mycobacterium_tuberculosis_h37rv', 'solanum_lycopersicum', 'neurospora_crassa', 'drosophila_melanogaster', 'geobacter_sulfurreducens_pca', 'sulfolobus_solfataricus_p2', 'cenarchaeum_symbiosum_a', 'streptococcus_pneumoniae_tigr4', 'gasterosteus_aculeatus', 'gallus_gallus', 'chlamydia_trachomatis_d_uw_3_cx', 'apis_mellifera', 'hyperthermus_butylicus_dsm_5456', 'borrelia_burgdorferi_b31', 'legionella_pneumophila_str_paris', 'homo_sapiens', 'pan_troglodytes', 'schizosaccharomyces_pombe', 'arabidopsis_thaliana', 'pyrobaculum_aerophilum_str_im2', 'caenorhabditis_elegans', 'archaeoglobus_fulgidus_dsm_4304', 'thermus_thermophilus_hb8', 'oryza_sativa', 'myxococcus_xanthus_dk_1622', 'methanobrevibacter_smithii_atcc_35061', 'aquifex_aeolicus_vf5', 'mycoplasma_pneumoniae_m129', 'thermoplasma_acidophilum_dsm_1728', 'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491', 'danio_rerio', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'enterococcus_faecalis_v583', 'ornithorhynchus_anatinus', 'halobacterium_salinarum_r1', 'physcomitrella_patens', 'haemophilus_influenzae_rd_kw20', 'candidatus_korarchaeum_cryptofilum_opf8', 'methanosarcina_acetivorans_c2a', 'yersinia_pestis_biovar_microtus_str_91001', 'staphylococcus_aureus_subsp_aureus_n315', 'brucella_abortus_bv_1_str_9_941', 'anolis_carolinensis', 'monodelphis_domestica', 'chloroflexus_aurantiacus_j_10_fl', 'aspergillus_nidulans', 'nanoarchaeum_equitans_kin4_m']
	dico = {'list': list,
	'Arc': ['methanosarcina_acetivorans_c2a', 'archaeoglobus_fulgidus_dsm_4304',
			'thermoplasma_acidophilum_dsm_1728', 'halobacterium_salinarum_r1',
			'methanobrevibacter_smithii_atcc_35061', 'pyrococcus_horikoshii_ot3',
			'nanoarchaeum_equitans_kin4_m', 'pyrobaculum_aerophilum_str_im2',
			'hyperthermus_butylicus_dsm_5456', 'candidatus_korarchaeum_cryptofilum_opf8',
			'cenarchaeum_symbiosum_a'],
	'ArcBac': ['geobacter_sulfurreducens_pca', 'pyrococcus_horikoshii_ot3',
			'yersinia_pestis_biovar_microtus_str_91001', 'thermus_thermophilus_hb8',
			'chloroflexus_aurantiacus_j_10_fl', 'halobacterium_salinarum_r1',
			'candidatus_korarchaeum_cryptofilum_opf8', 'archaeoglobus_fulgidus_dsm_4304',
			'myxococcus_xanthus_dk_1622', 'cenarchaeum_symbiosum_a',
			'brucella_abortus_bv_1_str_9_941', 'chlamydia_trachomatis_d_uw_3_cx',
			'methanosarcina_acetivorans_c2a', 'pyrobaculum_aerophilum_str_im2',
			'mycobacterium_tuberculosis_h37rv', 'neisseria_meningitidis_z2491',
			'hyperthermus_butylicus_dsm_5456', 'escherichia_coli_str_k_12_substr_mg1655',
			'legionella_pneumophila_str_paris', 'streptococcus_pneumoniae_tigr4',
			'aquifex_aeolicus_vf5', 'bacillus_subtilis_subsp_subtilis_str_168'],
	'ArcBacEuk': ['chlamydomonas_reinhardtii', 'solanum_lycopersicum', 'physcomitrella_patens',
			'vitis_vinifera', 'cenarchaeum_symbiosum_a', 'thermus_thermophilus_hb8',
			'archaeoglobus_fulgidus_dsm_4304', 'neurospora_crassa', 'aspergillus_nidulans',
			'myxococcus_xanthus_dk_1622', 'candidatus_korarchaeum_cryptofilum_opf8',
			'pongo_abelii', 'anolis_carolinensis', 'homo_sapiens', 'pan_troglodytes',
			'mus_musculus', 'gallus_gallus', 'monodelphis_domestica',
			'ornithorhynchus_anatinus', 'danio_rerio', 'methanosarcina_acetivorans_c2a',
			'gasterosteus_aculeatus', 'pyrococcus_horikoshii_ot3', 'caenorhabditis_elegans',
			'oryza_sativa', 'drosophila_melanogaster', 'arabidopsis_thaliana',
			'legionella_pneumophila_str_paris', 'geobacter_sulfurreducens_pca',
			'mycobacterium_tuberculosis_h37rv', 'brucella_abortus_bv_1_str_9_941',
			'apis_mellifera', 'leishmania_major', 'halobacterium_salinarum_r1',
			'pyrobaculum_aerophilum_str_im2', 'mycoplasma_pneumoniae_m129',
			'yersinia_pestis_biovar_microtus_str_91001', 'chloroflexus_aurantiacus_j_10_fl',
			'aquifex_aeolicus_vf5', 'anaplasma_phagocytophilum_str_hz',
			'neisseria_meningitidis_z2491', 'bacillus_subtilis_subsp_subtilis_str_168',
			'hyperthermus_butylicus_dsm_5456', 'streptococcus_pneumoniae_tigr4',
			'escherichia_coli_str_k_12_substr_mg1655', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819',
			'enterococcus_faecalis_v583', 'haemophilus_influenzae_rd_kw20'],
	'ArcEuk': ['chlamydomonas_reinhardtii', 'pongo_abelii', 'pan_troglodytes', 'mus_musculus',
			'anolis_carolinensis', 'ornithorhynchus_anatinus', 'monodelphis_domestica',
			'gallus_gallus', 'vitis_vinifera', 'homo_sapiens', 'danio_rerio', 'oryza_sativa',
			'caenorhabditis_elegans', 'candidatus_korarchaeum_cryptofilum_opf8',
			'physcomitrella_patens', 'pyrobaculum_aerophilum_str_im2', 'aspergillus_nidulans',
			'neurospora_crassa', 'gasterosteus_aculeatus', 'halobacterium_salinarum_r1',
			'archaeoglobus_fulgidus_dsm_4304', 'arabidopsis_thaliana', 'pyrococcus_horikoshii_ot3',
			'cenarchaeum_symbiosum_a', 'solanum_lycopersicum', 'drosophila_melanogaster',
			'methanosarcina_acetivorans_c2a', 'hyperthermus_butylicus_dsm_5456',
			'nanoarchaeum_equitans_kin4_m', 'leishmania_major', 'thermoplasma_acidophilum_dsm_1728',
			'apis_mellifera', 'dictyostelium_discoideum', 'methanobrevibacter_smithii_atcc_35061',
			'saccharomyces_cerevisiae', 'schizosaccharomyces_pombe'],
	'Bac': ['geobacter_sulfurreducens_pca', 'chloroflexus_aurantiacus_j_10_fl',
			'myxococcus_xanthus_dk_1622', 'thermus_thermophilus_hb8',
			'escherichia_coli_str_k_12_substr_mg1655', 'legionella_pneumophila_str_paris',
			'yersinia_pestis_biovar_microtus_str_91001', 'anaplasma_phagocytophilum_str_hz',
			'aquifex_aeolicus_vf5', 'chlamydia_trachomatis_d_uw_3_cx', 'bacillus_subtilis_subsp_subtilis_str_168',
			'brucella_abortus_bv_1_str_9_941', 'mycobacterium_tuberculosis_h37rv'],
	'BacEuk': ['mus_musculus', 'thermus_thermophilus_hb8', 'gallus_gallus',
			'geobacter_sulfurreducens_pca', 'chlamydomonas_reinhardtii',
			'monodelphis_domestica', 'danio_rerio', 'anolis_carolinensis',
			'pan_troglodytes', 'pongo_abelii', 'homo_sapiens', 'gasterosteus_aculeatus',
			'neurospora_crassa', 'leishmania_major', 'ornithorhynchus_anatinus',
			'solanum_lycopersicum', 'vitis_vinifera', 'oryza_sativa',
			'physcomitrella_patens', 'myxococcus_xanthus_dk_1622',
			'mycobacterium_tuberculosis_h37rv', 'anaplasma_phagocytophilum_str_hz',
			'caenorhabditis_elegans', 'drosophila_melanogaster', 'arabidopsis_thaliana',
			'apis_mellifera', 'chloroflexus_aurantiacus_j_10_fl', 'aquifex_aeolicus_vf5',
			'aspergillus_nidulans', 'legionella_pneumophila_str_paris',
			'mycoplasma_pneumoniae_m129', 'yersinia_pestis_biovar_microtus_str_91001',
			'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491',
			'chlamydia_trachomatis_d_uw_3_cx', 'bacillus_subtilis_subsp_subtilis_str_168',
			'dictyostelium_discoideum', 'brucella_abortus_bv_1_str_9_941'],
	'Euk': ['ornithorhynchus_anatinus', 'chlamydomonas_reinhardtii', 'mus_musculus',
			'vitis_vinifera', 'pongo_abelii', 'pan_troglodytes', 'homo_sapiens',
			'solanum_lycopersicum', 'gallus_gallus', 'anolis_carolinensis',
			'monodelphis_domestica', 'gasterosteus_aculeatus', 'physcomitrella_patens',
			'neurospora_crassa', 'arabidopsis_thaliana', 'danio_rerio', 'aspergillus_nidulans',
			'oryza_sativa', 'apis_mellifera', 'caenorhabditis_elegans', 'drosophila_melanogaster',
			'leishmania_major', 'schizosaccharomyces_pombe', 'dictyostelium_discoideum',
			'saccharomyces_cerevisiae']}
	dicoOrtho = {'list': list,
	'Arc': ['halobacterium_salinarum_r1', 'hyperthermus_butylicus_dsm_5456', 'methanobrevibacter_smithii_atcc_35061', 'candidatus_korarchaeum_cryptofilum_opf8', 'methanosarcina_acetivorans_c2a', 'pyrococcus_horikoshii_ot3', 'cenarchaeum_symbiosum_a', 'thermoplasma_acidophilum_dsm_1728', 'sulfolobus_solfataricus_p2', 'pyrobaculum_aerophilum_str_im2', 'nanoarchaeum_equitans_kin4_m', 'archaeoglobus_fulgidus_dsm_4304'],
	'ArcBac': ['bacillus_subtilis_subsp_subtilis_str_168', 'francisella_tularensis_subsp_tularensis_schu_s4', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'anaplasma_phagocytophilum_str_hz', 'pyrococcus_horikoshii_ot3', 'mycobacterium_tuberculosis_h37rv', 'geobacter_sulfurreducens_pca', 'sulfolobus_solfataricus_p2', 'cenarchaeum_symbiosum_a', 'streptococcus_pneumoniae_tigr4', 'chlamydia_trachomatis_d_uw_3_cx', 'hyperthermus_butylicus_dsm_5456', 'borrelia_burgdorferi_b31', 'legionella_pneumophila_str_paris', 'pyrobaculum_aerophilum_str_im2', 'archaeoglobus_fulgidus_dsm_4304', 'thermus_thermophilus_hb8', 'myxococcus_xanthus_dk_1622', 'aquifex_aeolicus_vf5', 'methanobrevibacter_smithii_atcc_35061', 'mycoplasma_pneumoniae_m129', 'thermoplasma_acidophilum_dsm_1728', 'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'enterococcus_faecalis_v583', 'halobacterium_salinarum_r1', 'haemophilus_influenzae_rd_kw20', 'candidatus_korarchaeum_cryptofilum_opf8', 'methanosarcina_acetivorans_c2a', 'yersinia_pestis_biovar_microtus_str_91001', 'staphylococcus_aureus_subsp_aureus_n315', 'brucella_abortus_bv_1_str_9_941', 'chloroflexus_aurantiacus_j_10_fl', 'nanoarchaeum_equitans_kin4_m'],
	'ArcBacEuk': ['bacillus_subtilis_subsp_subtilis_str_168', 'francisella_tularensis_subsp_tularensis_schu_s4', 'leishmania_major', 'pongo_abelii', 'saccharomyces_cerevisiae', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'anaplasma_phagocytophilum_str_hz', 'pyrococcus_horikoshii_ot3', 'mus_musculus', 'mycobacterium_tuberculosis_h37rv', 'solanum_lycopersicum', 'neurospora_crassa', 'geobacter_sulfurreducens_pca', 'sulfolobus_solfataricus_p2', 'cenarchaeum_symbiosum_a', 'streptococcus_pneumoniae_tigr4', 'gasterosteus_aculeatus', 'chlamydia_trachomatis_d_uw_3_cx', 'hyperthermus_butylicus_dsm_5456', 'borrelia_burgdorferi_b31', 'legionella_pneumophila_str_paris', 'homo_sapiens', 'pan_troglodytes', 'arabidopsis_thaliana', 'pyrobaculum_aerophilum_str_im2', 'nanoarchaeum_equitans_kin4_m', 'caenorhabditis_elegans', 'archaeoglobus_fulgidus_dsm_4304', 'thermus_thermophilus_hb8', 'oryza_sativa', 'myxococcus_xanthus_dk_1622', 'methanobrevibacter_smithii_atcc_35061', 'aquifex_aeolicus_vf5', 'mycoplasma_pneumoniae_m129', 'thermoplasma_acidophilum_dsm_1728', 'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491', 'danio_rerio', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'enterococcus_faecalis_v583', 'ornithorhynchus_anatinus', 'halobacterium_salinarum_r1', 'haemophilus_influenzae_rd_kw20', 'candidatus_korarchaeum_cryptofilum_opf8', 'methanosarcina_acetivorans_c2a', 'yersinia_pestis_biovar_microtus_str_91001', 'staphylococcus_aureus_subsp_aureus_n315', 'brucella_abortus_bv_1_str_9_941', 'anolis_carolinensis', 'monodelphis_domestica', 'chloroflexus_aurantiacus_j_10_fl', 'aspergillus_nidulans', 'gallus_gallus'],
	'Bac': ['bacillus_subtilis_subsp_subtilis_str_168', 'francisella_tularensis_subsp_tularensis_schu_s4', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'anaplasma_phagocytophilum_str_hz', 'mycobacterium_tuberculosis_h37rv', 'geobacter_sulfurreducens_pca', 'streptococcus_pneumoniae_tigr4', 'chlamydia_trachomatis_d_uw_3_cx', 'borrelia_burgdorferi_b31', 'legionella_pneumophila_str_paris', 'thermus_thermophilus_hb8', 'myxococcus_xanthus_dk_1622', 'aquifex_aeolicus_vf5', 'mycoplasma_pneumoniae_m129', 'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'enterococcus_faecalis_v583', 'haemophilus_influenzae_rd_kw20', 'yersinia_pestis_biovar_microtus_str_91001', 'staphylococcus_aureus_subsp_aureus_n315', 'brucella_abortus_bv_1_str_9_941', 'chloroflexus_aurantiacus_j_10_fl'],
	'BacEuk': ['bacillus_subtilis_subsp_subtilis_str_168', 'francisella_tularensis_subsp_tularensis_schu_s4', 'pongo_abelii', 'saccharomyces_cerevisiae', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'anaplasma_phagocytophilum_str_hz', 'mus_musculus', 'mycobacterium_tuberculosis_h37rv', 'neurospora_crassa', 'drosophila_melanogaster', 'geobacter_sulfurreducens_pca', 'streptococcus_pneumoniae_tigr4', 'gasterosteus_aculeatus', 'chlamydia_trachomatis_d_uw_3_cx', 'borrelia_burgdorferi_b31', 'legionella_pneumophila_str_paris', 'homo_sapiens', 'schizosaccharomyces_pombe', 'pan_troglodytes', 'arabidopsis_thaliana', 'caenorhabditis_elegans', 'thermus_thermophilus_hb8', 'oryza_sativa', 'myxococcus_xanthus_dk_1622', 'aquifex_aeolicus_vf5', 'mycoplasma_pneumoniae_m129', 'escherichia_coli_str_k_12_substr_mg1655', 'neisseria_meningitidis_z2491', 'danio_rerio', 'ornithorhynchus_anatinus', 'enterococcus_faecalis_v583', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'haemophilus_influenzae_rd_kw20', 'yersinia_pestis_biovar_microtus_str_91001', 'anolis_carolinensis', 'monodelphis_domestica', 'staphylococcus_aureus_subsp_aureus_n315', 'brucella_abortus_bv_1_str_9_941', 'aspergillus_nidulans', 'gallus_gallus'],
	'Euk': ['pongo_abelii', 'leishmania_major', 'saccharomyces_cerevisiae', 'chlamydomonas_reinhardtii', 'vitis_vinifera', 'mus_musculus', 'neurospora_crassa', 'solanum_lycopersicum', 'drosophila_melanogaster', 'gasterosteus_aculeatus', 'apis_mellifera', 'homo_sapiens', 'pan_troglodytes', 'schizosaccharomyces_pombe', 'arabidopsis_thaliana', 'caenorhabditis_elegans', 'oryza_sativa', 'danio_rerio', 'ornithorhynchus_anatinus', 'physcomitrella_patens', 'monodelphis_domestica', 'anolis_carolinensis', 'aspergillus_nidulans', 'gallus_gallus']}
	return dicoOrtho

def getDicoSpClean():
	dicoSp = {'anaplasma_phagocytophilum_str_hz': 'Anaplasmaphagocytophilum',
	'anolis_carolinensis': 'Anoliscarolinensis',
	'apis_mellifera': 'Apismellifera',
	'aquifex_aeolicus_vf5': 'Aquifexaeolicus',
	'arabidopsis_thaliana': 'Arabidopsisthaliana',
	'archaeoglobus_fulgidus_dsm_4304': 'Archeoglobusfulgidus',
	'aspergillus_nidulans': 'Aspergillusnidulans',
	'bacillus_subtilis_subsp_subtilis_str_168': 'Bacillussubtilis',
	'borrelia_burgdorferi': 'Borreliaburgdorferi',
	'brucella_abortus_bv_1_str_9_941': 'Brucellaabortus',
	'caenorhabditis_elegans': 'Caenorhabditiselegans',
	'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819': 'Campylobacterjejuni',
	'candidatus_korarchaeum_cryptofilum_opf8': 'Candidatuskorarchaeumcryptofilum',
	'cenarchaeum_symbiosum_a': 'Cenarchaeumsymbiosum',
	'chlamydia_trachomatis_d_uw_3_cx': 'Chlamydiatrachomatis',
	'chlamydomonas_reinhardtii': 'Chlamydomonasreinhardtii',
	'chloroflexus_aurantiacus_j_10_fl': 'Chloroflexusaurantiacus',
	'danio_rerio': 'Daniorerio',
	'dictyostelium_discoideum': 'Dictyosteliumdiscoideum',
	'drosophila_melanogaster': 'Drosophilamelanogaster',
	'enterococcus_faecalis_v583': 'Enterococcusfaecalis',
	'escherichia_coli_str_k_12_substr_mg1655': 'Escherichiacoli',
	'francisella_tularensis_subsp_tularensis_schu_s4': 'Francisellatularensis',
	'gallus_gallus': 'Gallusgallus',
	'gasterosteus_aculeatus': 'Gasterosteusaculeatus',
	'geobacter_sulfurreducens_pca': 'Geobactersulfurreducens',
	'haemophilus_influenzae_rd_kw20': 'Haemophilusinfluenzae',
	'halobacterium_salinarum_r1': 'Halobacteriumsalinarum',
	'homo_sapiens': 'Homosapiens',
	'hyperthermus_butylicus_dsm_5456': 'Hyperthermusbutylicus',
	'legionella_pneumophila_str_paris': 'Legionellapneumoniae',
	'leishmania_major': 'Leishmaniamajor',
	'methanobrevibacter_smithii_atcc_35061': 'Metanobrevibactersmithii',
	'methanosarcina_acetivorans_c2a': 'Methanosarcinaacetivorans',
	'monodelphis_domestica': 'Monodelphisdomestica',
	'mus_musculus': 'Musmusculus',
	'mycobacterium_tuberculosis_h37rv': 'Mycobacteriumtuberculosis',
	'mycoplasma_pneumoniae_m129': 'Mycoplasmapneumoniae',
	'myxococcus_xanthus_dk_1622': 'Myxococcusxanthus',
	'nanoarchaeum_equitans_kin4_m': 'Nanoarcheumequitans',
	'neisseria_meningitidis_z2491': 'Neisseriameningitidis',
	'neurospora_crassa': 'Neurosporacrassa',
	'ornithorhynchus_anatinus': 'Ornithorhynchusanatinus',
	'oryza_sativa': 'Oryzasativa',
	'pan_troglodytes': 'Pantroglodytes',
	'physcomitrella_patens': 'Physcomitrellapatens',
	'pongo_abelii': 'Pongoabelii',
	'pyrobaculum_aerophilum_str_im2': 'Pyrobaculumaerophilum',
	'pyrococcus_horikoshii_ot3': 'Pyrococcushorikoshii',
	'saccharomyces_cerevisiae': 'Saccharomycescerevisiae',
	'schizosaccharomyces_pombe': 'Schizosaccharomycespombe',
	'solanum_lycopersicum': 'Solanumlycopersicum',
	'staphylococcus_aureus_subsp_aureus_n315': 'Staphylococcusaureus',
	'streptococcus_pneumoniae_tigr4': 'Streptococcuspneumoniae',
	'sulfolobus_solfataricus': 'Sulfolobussolfataricus',
	'thermoplasma_acidophilum_dsm_1728': 'Thermoplasmaacidophilum',
	'thermus_thermophilus_hb8': 'Thermusthermophilus',
	'vibrio_cholerae': 'Vibriocholerae',
	'vitis_vinifera': 'Vitisvinifera',
	'wolbachia_endosymbiont_of_drosophila_melanogaster': 'Wolbachiaendosymbiontofdrosophilamelanogaster',
	'yersinia_pestis_biovar_microtus_str_91001': 'Yersiniapestis',
	'sulfolobus_solfataricus_p2': 'Sulfolobussolfataricus',
	'borrelia_burgdorferi_b31': 'Borreliaburgdorferi'}
	return dicoSp

def importTreeDicoQuality():
	dico = {'Euk': ['8005', '19494', '1542', '3765', '9449', '8110', '7899', '3064', '1023', '9660', '7021', '9214', '14798', '6309', '20015', '17206', '19644', '5879', '15878', '4280', '3018', '4570', '2070', '5601', '3729', '5183', '3056', '8241', '5481', '3900', '15946', '8136', '14829', '11076', '1805', '16772', '6475', '6876', '6816', '7187', '11232', '9270', '5563', '1454', '2139', '5778', '16731', '6953', '1273', '4432', '9458', '7250', '14985', '4723', '18476', '2041', '15330', '1705', '7261', '5723', '9195', '20328', '15378', '4364', '8135', '5242', '16631', '19311', '8861', '1792', '746', '2885', '4415', '9032', '9325', '3046', '4195', '3909', '9682', '5254'],
    'Bac': ['10039', '10124', '13877', '14220', '398', '10147', '13950', '411', '312', '202', '9819', '10056', '13949', '10197', '7886', '7940', '10031', '14030', '458', '13910', '13966', '302', '9980', '14456', '10192', '311', '13925', '342', '38', '9964', '7991', '13', '10130', '171', '9927', '113', '228', '9919', '9817', '13976', '205', '7982', '9900', '10144', '9914', '10182', '379', '10167', '15525', '8013', '29', '284', '56', '191', '439', '10030', '344', '7883', '14255', '14071', '438', '89', '469', '7957', '206', '264', '9882', '10161', '14544', '281', '14219', '14226', '13944', '10159', '212', '13936', '14365', '14395', '352', '10150', '288', '7861', '88', '9897', '486', '14281', '9', '14362', '392', '20954', '14346', '14457', '13823', '14260', '14223', '14057', '9886', '15558', '14119', '14123', '13885', '14699', '9951', '9888', '424', '14259', '9996', '7889', '10022', '145', '7916', '357', '13824', '14590', '10143', '7984', '7917', '10177', '9949', '13799', '14673', '10114', '10113', '10033', '14153', '9829', '9958', '13929', '9976', '9893', '10141', '10', '14313', '305', '14060', '8016', '220', '20990'],
    'Arc': ['13473', '13351', '13118', '13235', '13322', '13313', '13420', '13521', '13355', '13526', '13174', '13292', '13306', '13529', '15775', '13547', '13275', '13341', '13429', '13300', '26714', '13516', '15654', '13384', '13162', '13430', '13111', '13354', '13180', '13387', '13271', '15708', '13145', '13119', '13129', '13391', '13237', '13426', '15699', '13261', '13183', '13340', '13273', '15714', '13472', '13260', '13187', '13325', '13419', '15726', '13182', '13392', '13487', '13522', '13184', '13139', '13267', '13152', '13476', '13385', '13502', '13368', '13530', '13359', '13236', '13229', '13323', '15715', '13495', '13311', '13120', '13501', '13208', '13193', '13439'],
	'ArcBac': ['10069', '453', '13838', '14050', '14169', '10217', '14603', '167', '139', '111', '354', '13251', '307', '10051', '13246', '14082', '14173', '10176', '13356', '10040', '13906', '9950', '9977', '10227', '20956', '7885', '484', '27', '10027', '68', '9884', '9898', '10136', '230', '100', '14297', '408', '9895', '10149', '461', '74', '14487', '293', '9936', '9867', '72', '7955', '10129', '15784', '14643', '13245', '9930', '10034', '227', '13423', '13498', '413', '13232', '10198', '14448', '13829', '9828', '10068', '10172', '13893', '14492', '10001', '10209', '10232', '13343', '7911', '13546', '9929', '13175', '13332', '214', '15710', '9844', '9857', '9955', '10041', '372', '10194', '14214', '14245', '9969', '13134', '13176', '14231', '7962', '13216', '456', '14446', '9833', '13404', '13144', '13447', '13458', '10066', '8053', '13422', '14035', '10230', '9851', '13952', '10042', '489', '7987', '14193', '10013', '14085', '28', '9913', '14618', '10058', '14191', '14534', '9982', '13542', '13541', '13981', '13543', '13480', '13367', '10231', '10156', '13192', '15785', '21178', '13223', '497', '380', '14136', '8040', '10166', '10050', '482', '15782', '10128', '14177', '7882', '13828', '14484', '13179', '14536', '301', '13485', '10118', '9904', '10223', '13441', '14383', '14217', '13257', '7902', '62', '16535', '21271', '9869', '13248', '14324', '9834', '15709', '13244', '14109', '13171', '10016', '14162', '494', '13167', '339', '9879', '9881', '10226', '140', '13225', '84', '10102', '13886', '13262', '14248', '13151', '9874', '9988', '332', '13407', '14449', '427', '14599', '7875', '14187', '13438', '13326', '14012', '14520', '15758', '266', '13282', '14662', '13269', '13433', '13297', '10059', '13403', '10174', '9915', '244', '13933', '15630', '94', '13258', '13159', '10215', '10052', '9905', '10169', '7853', '13298'],
	'ArcBacEuk': ['502', '196', '198', '300', '369', '225', '13115', '10196', '329', '10153'],
	'BacEuk': ['10076', '11', '136']}
	return dico

def main(path):
	listSp = getSpParse()
	geniListDir = path + 'Homology/'
	dicoTreeGododQual = importTreeDicoQuality()
	dicoGeSp = importGeneList(geniListDir, listSp['list'])

	grpSp = {'Arc': [], 'Bac': [], 'Euk': [], 'ArcBac': [], 'ArcBacEuk': [], 'BacEuk': []}

	dicoDF = {}
	indexList = []
	for sp in listSp['list']:
		dicoDF[sp]=[]
	cptNbFam = 0
	for grp in grpSp:
		print(grp)
		for tree in dicoTreeGododQual[grp]:
			if grp == 'BacEuk':
				inputpG4Family = path+"mergedpG4FamiliesPecentrOrtho/"+grp+'/Families_'+tree+".csv"
			else:
				inputpG4Family = path+"mergedpG4FamiliesPecentrOrtho/"+grp+'/'+tree+".csv"
			try:
				with open(inputpG4Family) as f:
					content = f.read()
			except:
				print('could not open '+inputpG4Family)
			else:
				with open(inputpG4Family) as f:
					content = f.read()
					lines = content.split('\n')
					lines = lines[1:]
					if lines[-1] == '':
						lines = lines[:-1]
					nbFam = len(lines)
					for l in lines:
						for sp in dicoDF:
							dicoDF[sp].append(0)
						w = l.split('\t')
						fam = w[0]
						if fam != '':
							i = tree+'_'+str(cptNbFam)
							indexList.append(i)
							genes = w[3].split('|')
							if genes[-1] == '':
								genes = genes[:-1]
							for g in genes:
								spFam = dicoGeSp[g]
								dicoDF[spFam][cptNbFam] += 1
							cptNbFam += 1
						else:
							for sp in dicoDF:
								dicoDF[sp] = dicoDF[sp][:-1]
							# indexList = indexList[:-1]
		dicoNameSp = getDicoSpClean()
		dicoDfFin = {}
		for sp in dicoDF:
			dicoDfFin[dicoNameSp[sp]] = dicoDF[sp]
		df = pd.DataFrame(data=dicoDfFin)
		df.index = indexList
		df.to_csv(path_or_buf=path+'Matrices/AllTreeOrtho_Extended'+grp+'.csv', header=True, index=True, sep='\t')

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getMatricepG4Fan')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	main(path)

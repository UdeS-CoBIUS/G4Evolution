---
title: "Readme"
author: "Anaïs Vannutelli"
output: html_document
---

## Prerequesit

Librairy used:

* Snakemake
* Python 3.9.5
* pandas 1.2.4
* seaborn 0.11.1
* matplotlib 3.4.2
* numpy 1.20.3
* jupyter notebook 6.4.0
* biopython 1.78
* networkx 2.5.1

## Main scripts description

* **orthologyGraph.py** &rarr; from EnsemblCompara homology file (https://ftp.ensemblgenomes.ebi.ac.uk/pub/pan_ensembl/release-46/tsv/ensembl-compara/homologies/Compara.99.protein_default.homologies.tsv.gz) create a graph to retrieve groups of orthologs 1 to 1. This script creates one file per group, and output in it genes name from the group. The name group is randomly provided.
* **getOrthoFasta.py** &rarr; from the files containing 1:1 orthologous genes, create a fasta file containing all gene sequences from these groups of orthologs.
* **launchkAlignMultiprocess.py** &rarr; using the fasta files of 1:1 orthologous genes, launch kAlign to get multiple sequence alignments.
* **processAlignment.py** &rarr; this script aims to read one MSA and parsed it into a matrice corresponding to pG4: pG4 correspond to 0, pG4 correspond to 1 and gaps correspond to -.
* **getpG4Family.py** &rarr; from the previously generated matrice, the matrice is browsed to get overlaping pG4 by retrieving segments of the alignments with ones. One segment correspond to one pG4 family.
* **runQGRSMapperOnFam.py** &rarr; for eahc pG4 family, retrieve the corresponding alignment sequences to run QGRSmapper and thus predict G4 by homology. G4RNA screener is also used to get which G4 were originally predicted.
* **generateShuffle.py** &rarr; from a list of orthologous genes (those with a good alignments), generate 10 mono-nucleotids shuffled sequences, 10 di-nucleotids shuffled sequences and 10 tri-nucleotids shuffled sequences.
* **G4AnnoationShuffle.py** &rarr; from G4RNAscreener ouput, find out where are the pG4 in the shuffle dataset.
* **processAlignmentShuffle.py**, **getpG4FamilyShuffle.py** and **runQGRSMapperOnFamShuffle.py** are the same scripts than the non shuffle ones. Some modifications have been made to make them run with shuffled sequences, in particular due the identification names which are different.

## Jupyter notebook scripts:

* **ShuffleRepro.ipynb** &rarr; analyse results from the shuffle dataset.
* **OrthologyGrpAnalyse.ipynb** &rarr; analyse orthology groups and alignments.

## Usefull scripts:

* **convertXMLtoNewick.py** &rarr; convert an XML phylogenetic tree to a newick format.
* **getAlignmentPosi.py** &rarr; take as input a filename of an alignment, a start and an end location, the output is the alginment sequences corresponding to these information.
* **getBigMatrices.py** &rarr; to use Count, matrices are requiered to known in which species are present/absent pG4 families. This script is creating a matrice for each species clades.


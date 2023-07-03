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

## Script description

* **BAM_correspondance_Prot-Repl.py** &rarr; parse the tsv file from ENCORE to download all RBPs eCLIP data.
* **generateG4RandomCLIP.py** &rarr; generate a fake pG4 dataset from predicted G4 contained in GAIA db from the first release in June 2022. This script start to get the number n of transcripts, then a list of n 0 is created and for the number of pG4, a random element of the list is implemented by a +1. Then, for each random pG4, a random location between the start and end of the transcript is selected and allocated to the pG4. Thus the number of pG4 is the same between the normal dataset and the random one, but their location is different.
* **ResulteCLIP.ipynb** &rarr; jupyter notebook script to analyse results and generate figures.

## Snakefile workflow

Two snakemake workflow were made: one for humans and on for POSTAR species (mouse, fly, S. cerevisiae, zebra fish and worm). These workflow are making the same important steps which is using bedtool to compare bed files of pG4 and RNA binding sites of RBPs.

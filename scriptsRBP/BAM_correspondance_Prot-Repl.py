# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd

def main(path):
    with open(path+'data/eCLIP/bam/urlFile.txt') as f:
        listBamFile = f.read().split('\n')
    del listBamFile[0]
    del listBamFile[-1]
    listBamFile = [x.split('/')[4] for x in listBamFile if '.bam' in x]

    with open(path+'data/eCLIP/bam/urlFileControl.txt') as f:
        listBamFileCtrl = f.read().split('\n')
    del listBamFileCtrl[0]
    del listBamFileCtrl[-1]
    listBamFileCtrl = [x.split('/')[4] for x in listBamFileCtrl if '.bam' in x]

    dfExp = pd.read_csv(path+'data/eCLIP/bam/experiment_report_2022_5_6_16h_25m.tsv',
        sep='\t', skiprows=1)

    dfExpCtrl = pd.read_csv(path+'data/eCLIP/bam/experiment_report_CONTROL_2022_5_10_12h_13m.tsv',
        sep='\t', skiprows=1)

    dicoCorrespondance = {'Protein': [], 'cellType': [], 'Species': [], 'Experiment': [], 'BamName': []}

    for index, row in dfExp.iterrows():
        listFiles = row['Files'].split(',')
        BamFile = [x.split('/')[2] for x in listFiles if x.split('/')[2] in listBamFile]
        for f in BamFile:
            dicoCorrespondance['Protein'].append(row['Target gene symbol'])
            dicoCorrespondance['cellType'].append(row['Biosample term name'])
            dicoCorrespondance['Species'].append(row['Organism'])
            dicoCorrespondance['Experiment'].append(row['Assay title'])
            dicoCorrespondance['BamName'].append(f)

    for index, row in dfExpCtrl.iterrows():
        listFiles = row['Files'].split(',')
        BamFile = [x.split('/')[2] for x in listFiles if x.split('/')[2] in listBamFileCtrl]
        for f in BamFile:
            dicoCorrespondance['Protein'].append(row['Description'].split('against ')[-1])
            dicoCorrespondance['cellType'].append(row['Biosample term name'])
            dicoCorrespondance['Species'].append(row['Organism'])
            dicoCorrespondance['Experiment'].append(row['Assay title'])
            dicoCorrespondance['BamName'].append(f)

    dfCorrespondance = pd.DataFrame(data=dicoCorrespondance)
    dfCorrespondance.to_csv(path_or_buf=path + 'data/eCLIP/bam/FileExperimentBamName.csv', header=True, index=None, sep='\t')

    outputFile = []
    groups = dfCorrespondance.groupby(['Protein','cellType'])
    for name, group in groups:
        if len(list(group[group.Experiment == 'eCLIP']['BamName'])) > 0 and \
            'Control eCLIP' in list(group.Experiment):
            # print(group)
            cmd = "Rscript /home/vana2406/scratch/G4RBPevo/scripts/nearByndingWorkflow.R "+\
                str(list(group[group.Experiment == 'eCLIP']['BamName'])[0])+" "+\
                str(list(group[group.Experiment == 'eCLIP']['BamName'])[1])+" "+\
                name[1]+" "+\
                str(list(group[group.Experiment == 'Control eCLIP']['BamName'])[0])+" "+\
                name[0]
            outputFile.append(cmd)
    output = open('list_ghost_tasks_runNearBynding.txt', "w")
    output.write('\n'.join(outputFile))
    output.close()

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'getFastaTree')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)

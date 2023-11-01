import pandas as pd
import os
import itertools

###### Config file and sample sheets #####
# configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
samples['Treatment'] = [i[0] for i in samples.index.str.split('-')]
Treatments = samples['Treatment'].unique()
NonControlTreatments = [i for i in Treatments if i != 'DMSO']

autosomes = [str(i) for i in range(1,23)]

SAMPLE_NAMES = samples.index
AUTOSOMES = ['chr' + str(i) for i in range(1,23)]
CHROMS = AUTOSOMES + ['chrX', 'chrY']

titration_series_samples = pd.read_csv("config/samples.titrationseries.tsv", sep='\t', index_col=0)
chRNA_samples = pd.read_csv("config/samples.chRNAPilotSequencing.tsv", sep='\t', index_col=0)
ExpOf52_samples = pd.read_csv("config/samples.52MoleculeExperiment.tsv", sep='\t', index_col=0)
Exp_202310_3MoleculesOfInterest_samples = pd.read_csv("config/samples.3MoleculesOfInterestExperiment.tsv", sep='\t', index_col=0)

ExpOf52_samples_Treatments = ExpOf52_samples['Treatment'].unique()
ExpOf52_samples_NonControlTreatments = [i for i in ExpOf52_samples_Treatments if i != 'DMSO'] 

AllSamples = list(itertools.chain(*[i.index.unique().tolist() for i in [samples, titration_series_samples, chRNA_samples, ExpOf52_samples]]))
AllNEBNextSamples = list(itertools.chain(*[i.index.unique().tolist() for i in [titration_series_samples, chRNA_samples, ExpOf52_samples, Exp_202310_3MoleculesOfInterest_samples]]))
AllSamples_202310 = list(itertools.chain(*[i.index.unique().tolist() for i in [samples, titration_series_samples, chRNA_samples, ExpOf52_samples, Exp_202310_3MoleculesOfInterest_samples]]))


def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 52000

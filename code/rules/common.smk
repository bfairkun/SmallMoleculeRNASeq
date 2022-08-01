import pandas as pd
import os

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

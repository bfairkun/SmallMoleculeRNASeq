import pandas as pd
import os

###### Config file and sample sheets #####
# configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.tsv",sep='\t', index_col=0)
samples['Treatment'] = [i[0] for i in samples.index.str.split('-')]
Treatments = samples['Treatment'].unique()
NonControlTreatments = [i for i in Treatments if i != 'DMSO']

SAMPLE_NAMES = samples.index
AUTOSOMES = ['chr' + str(i) for i in range(1,23)]
CHROMS = AUTOSOMES + ['chrX', 'chrY']



    # TreatmentSamples = samples.loc[samples['Treatment']==treatment].index
    # for sample in TreatmentSamples:
    #     print(f"Something.{sample}.bam\t{treatment}")
    # for sample in ControlSamples:
    #     print(f"Something.{sample}.bam\tDMSO")

# # How to access values in samples.tsv

# print(samples)
# print( expand("Hello {sample}", sample=samples.index) )
# print( samples.at["A", "R1"] )

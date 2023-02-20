
rule GatherExpOf52Results:
    input:
        expand("FastqFastp/{sample}.R1.fastq.gz", sample=ExpOf52_samples.index),
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=ExpOf52_samples.index.tolist()),
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind.counts.gz",
        "featureCounts/AllSamples_Counts.txt",
        expand("SplicingAnalysis/leafcutter/groupsfiles/ExpOf52_{treatment}.txt", treatment=ExpOf52_samples_NonControlTreatments),
        "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab",
        expand("SplicingAnalysis/leafcutter/differential_splicing/ExpOf52_{treatment}_effect_sizes.txt", treatment=ExpOf52_samples_NonControlTreatments),

use rule CopyAndMergeFastq as CopyAndMergeFastq_ExpOf52 with:
    input:
        R1 = lambda wildcards: ExpOf52_samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: ExpOf52_samples.loc[wildcards.sample]['R2'],
    wildcard_constraints:
        sample = "|".join(ExpOf52_samples.index)


rule MakeGroupsFiles_ExpOf52:
    output:
        "SplicingAnalysis/leafcutter/groupsfiles/ExpOf52_{treatment}.txt"
    wildcard_constraints:
        treatment = "|".join(ExpOf52_samples_NonControlTreatments)
    run:
        with open(output[0], 'w') as fh:
            TreatmentSamples = ExpOf52_samples.loc[ExpOf52_samples['Treatment']==wildcards.treatment].index
            for sample in TreatmentSamples:
                _ = fh.write(f"{sample}\t{wildcards.treatment}\n")
            for sample in ExpOf52_samples.loc[ExpOf52_samples['Treatment']=='DMSO'].index:
                _ = fh.write(f"{sample}\tDMSO\n")

use rule leafcutter_ds as leafcutter_ds_ExpOf52 with:
    input:
        groupfile = "SplicingAnalysis/leafcutter/groupsfiles/ExpOf52_{treatment}.txt",
        numers = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind.counts.gz",
        Rscript = "scripts/leafcutter/scripts/leafcutter_ds.R"
    output:
        "SplicingAnalysis/leafcutter/differential_splicing/ExpOf52_{treatment}_effect_sizes.txt",
        "SplicingAnalysis/leafcutter/differential_splicing/ExpOf52_{treatment}_cluster_significance.txt"
    wildcard_constraints:
        treatment = "|".join(ExpOf52_samples_NonControlTreatments)
    params:
        Prefix = "SplicingAnalysis/leafcutter/differential_splicing/ExpOf52_",
        ExtraParams = "-i 2 -g 2"
    log:
        "logs/leafcutter_ds_ExpOf52/{treatment}.log"


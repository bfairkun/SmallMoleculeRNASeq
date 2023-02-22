
rule GatherExpOf52Results:
    input:
        expand("FastqFastp/{sample}.R1.fastq.gz", sample=ExpOf52_samples.index),
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=ExpOf52_samples.index.tolist()),
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind.counts.gz",
        "featureCounts/AllSamples_Counts.txt",
        expand("SplicingAnalysis/leafcutter/groupsfiles/ExpOf52_{treatment}.txt", treatment=ExpOf52_samples_NonControlTreatments),
        "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        expand("SplicingAnalysis/leafcutter/differential_splicing/ExpOf52_{treatment}_effect_sizes.txt", treatment=ExpOf52_samples_NonControlTreatments),
        "DE_testing/ExpOf52_Counts.mat.txt.gz",
        expand("bigwigs/unstranded/{sample}.bw", sample=AllSamples),
        expand("bigwigs/stranded/{sample}.minus.bw", sample=AllNEBNextSamples),
        expand("bigwigs/stranded/{sample}.plus.bw", sample=AllNEBNextSamples),
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numersJunctionCounts.bed.gz.tbi" ,
        "SplicingAnalysis/MergedExp52_Contrast/effect_sizes.txt",

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

rule MakeGroupsFiles_MergedExpOf52:
    output:
        "SplicingAnalysis/MergedExp52_Contrast/groups.txt"
    run:
        with open(output[0], 'w') as fh:
            for sample in ExpOf52_samples.loc[ExpOf52_samples['Treatment']!='DMSO'].index:
                _ = fh.write(f"{sample}\tTreated\n")
            for sample in ExpOf52_samples.loc[ExpOf52_samples['Treatment']=='DMSO'].index:
                _ = fh.write(f"{sample}\tDMSO\n")


use rule leafcutter_ds as leafcutter_ds_ExpOf52 with:
    input:
        groupfile = "SplicingAnalysis/leafcutter/groupsfiles/ExpOf52_{treatment}.txt",
        numers = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.counts.gz",
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

rule leafcutter_ds_ExpOf52_MergedContrast:
    input:
        groupfile = "SplicingAnalysis/MergedExp52_Contrast/groups.txt",
        numers = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.counts.gz",
        Rscript = "scripts/leafcutter/scripts/leafcutter_ds.R"
    output:
        "SplicingAnalysis/MergedExp52_Contrast/effect_sizes.txt",
        "SplicingAnalysis/MergedExp52_Contrast/cluster_significance.txt"
    params:
        Prefix = "SplicingAnalysis/MergedExp52_Contrast/",
        ExtraParams = ""
    log:
        "logs/leafcutter_ds_ExpOf52_MergedContrast.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript {input.Rscript} -p {threads} -o {params.Prefix} {params.ExtraParams} {input.numers} {input.groupfile} &> {log}
        """

rule DE_testing_ExpOf52:
    input:
        "featureCounts/AllSamples_Counts.txt"
    output:
        results = "DE_testing/ExpOf52_Results.txt.gz",
        counts = "DE_testing/ExpOf52_Counts.mat.txt.gz"
    log:
        "logs/DE_testing_ExpOf52.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/DE_edgeR_ExpOf52.R {input} {output.results} {output.counts} &> {log}
        """

rule MakePSITable:
    input:
        juncs = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.counts.gz"
    output:
        PSI = temp("SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.bed"),
        counts =temp("SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numersJunctionCounts.bed" )
    params:
        Prefix = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers"
    log:
        "logs/MakePSITable.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/leafcutter_to_PSI.R {input} {params.Prefix} &> {log}
        """

rule bgzip_and_tabix_bed:
    input:
        PSI = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.bed",
        counts ="SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numersJunctionCounts.bed" 
    output:
        PSI = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.bed.gz",
        counts ="SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numersJunctionCounts.bed.gz" ,
        PSI_tab = "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.bed.gz.tbi",
        counts_tab ="SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numersJunctionCounts.bed.gz.tbi" 
    shell:
        """
        bgzip {input.PSI}
        bgzip {input.counts}
        tabix -p bed {output.PSI}
        tabix -p bed {output.counts}
        """


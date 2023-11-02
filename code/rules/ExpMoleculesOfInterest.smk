
rule GatherExpMoleculesOfInterest:
    input:
        "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers.counts.gz",
        "featureCounts/AllSamples_Counts.txt",
        expand("SplicingAnalysis/Exp202310_3Molecules_Contrasts/{ContrastName}_cluster_significance.txt", ContrastName=Exp_202310_3MoleculesOfInterest_contrasts.index)

rule CopyAndMergeFastq_Exp_3MoleculesOfInterest:
    """
    Useful for when a single sample is spread across multiple fastq
    """
    input:
        R1 = lambda wildcards: Exp_202310_3MoleculesOfInterest_samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: Exp_202310_3MoleculesOfInterest_samples.loc[wildcards.sample]['R2'],
    output:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    wildcard_constraints:
        sample = "|".join(Exp_202310_3MoleculesOfInterest_samples.index)
    shell:
        """
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        """

use rule make_leafcutter_juncfile as make_leafcutter_juncfile_all_202310 with:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=AllSamples_202310),
    output:
        "SplicingAnalysis/leafcutter_all_samples_202310/juncfilelist.txt"

use rule leafcutter_cluster as leafcutter_cluster_all_202310 with:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=AllSamples_202310),
        juncfile_list = "SplicingAnalysis/leafcutter_all_samples_202310/juncfilelist.txt"
    output:
        "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers.counts.gz"
    params:
        "-r SplicingAnalysis/leafcutter_all_samples_202310/"
    log:
        "logs/leafcutter_cluster_all_202310.log"

use rule MakePSITable as MakePSITable_202310 with:
    input:
        juncs = "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers.counts.gz"
    output:
        PSI = temp("SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers.bed"),
        counts =temp("SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numersJunctionCounts.bed" )
    params:
        Prefix = "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers"
    log:
        "logs/MakePSITable_202310.log"

rule MakeGroupsFiles_ExpOf3:
    output:
        expand("SplicingAnalysis/Exp202310_3Molecules_Contrasts/{ContrastName}.groups.tsv", ContrastName = Exp_202310_3MoleculesOfInterest_contrasts.index)
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/Exp202310_MakeLeafcutterContrasts.R
        """

use rule leafcutter_ds as leafcutter_ds_ExpOf3 with:
    input:
        groupfile = "SplicingAnalysis/Exp202310_3Molecules_Contrasts/{treatment}.groups.tsv",
        numers = "SplicingAnalysis/leafcutter_all_samples_202310/leafcutter_perind_numers.counts.gz",
        Rscript = "scripts/leafcutter/scripts/leafcutter_ds.R"
    output:
        "SplicingAnalysis/Exp202310_3Molecules_Contrasts/{treatment}_effect_sizes.txt",
        "SplicingAnalysis/Exp202310_3Molecules_Contrasts/{treatment}_cluster_significance.txt"
    wildcard_constraints:
        ContrastName = "|".join(Exp_202310_3MoleculesOfInterest_contrasts.index)
    params:
        Prefix = "SplicingAnalysis/Exp202310_3Molecules_Contrasts/",
        ExtraParams = "-i 2 -g 2"
    log:
        "logs/leafcutter_ds_ExpOf3/{treatment}.log"

# rule TidyDoseResponseData_ExpOf3:
#     input:

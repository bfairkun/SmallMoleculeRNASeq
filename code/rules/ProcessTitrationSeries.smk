rule GatherTitrationResults:
    input:
        expand("Fastq/{sample}.R1.fastq.gz", sample=titration_series_samples.index),
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=titration_series_samples.index)

rule CopyAndMergeFastq:
    """
    Useful for when a single sample is spread across multiple fastq
    """
    input:
        R1 = lambda wildcards: titration_series_samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: titration_series_samples.loc[wildcards.sample]['R2'],
    output:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    wildcard_constraints:
        sample = "|".join(titration_series_samples.index)
    shell:
        """
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        """


use rule fastp as fastp_input_fn_from_wildcards with:
    input:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    wildcard_constraints:
        sample = "|".join(pd.concat([titration_series_samples, chRNA_samples, ExpOf52_samples, Exp_202310_3MoleculesOfInterest_samples]).index)


use rule featurecounts as featurecounts_titrationseries with:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= titration_series_samples.index.unique()),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai", sample= titration_series_samples.index.unique()),
        Comprehensive_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    log:
        "logs/featureCounts_titration_series.log"
    output:
        "featureCounts/Counts.titration_series.txt"
    params:
        extra = "-s 2 -p"

use rule make_leafcutter_juncfile as make_leafcutter_juncfile_titrationseries with:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=titration_series_samples.index.unique()),
    output:
        "SplicingAnalysis/leafcutter/juncfilelist.autosomes.titrationseries.txt"

use rule leafcutter_cluster as leafcutter_cluster_titrationseries with:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=titration_series_samples.index.unique()),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.titrationseries.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering/autosomes_titrationseries/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering/autosomes_titrationseries/leafcutter_perind_numers.counts.gz"
    log:
        "logs/leafcutter_cluster_titrationseries/autosomes.log"
    params:
        "-r SplicingAnalysis/leafcutter/clustering/autosomes_titrationseries/"

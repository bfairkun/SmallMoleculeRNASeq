
rule GatherChRNAResults:
    input:
        expand("FastqFastp/{sample}.R1.fastq.gz", sample=chRNA_samples.index),
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=chRNA_samples.index),
        expand("FragLenths/{sample}.1M.test.txt.gz", sample=pd.concat([titration_series_samples, chRNA_samples]).index),
        "Multiqc/multiqc_report.html",
        "featureCounts/Counts.chRNA.txt",
        "../output/QC/ReadCountsAndJunctionsPerSamples.tsv",
        expand("QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt", sample=chRNA_samples.index),
        "featureCounts/Counts.chRNA.txt",
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.counts.gz"

use rule CopyAndMergeFastq as CopyAndMergeFastq_chRNA with:
    input:
        R1 = lambda wildcards: chRNA_samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: chRNA_samples.loc[wildcards.sample]['R2'],
    output:
        R1 = temp("Fastq/{sample}.R1.fastq.gz"),
        R2 = temp("Fastq/{sample}.R2.fastq.gz"),
    wildcard_constraints:
        sample = "|".join(chRNA_samples.index)

rule CountReadsPerSample:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=AllSamples),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=AllSamples),
    output:
        "../output/QC/ReadCountsAndJunctionsPerSamples.tsv"
    log:
        "logs/CountReadsPerSample.log"
    shell:
        """
        # exec > {log} 2>&1
        # set -x
        for f in {input.bam}
        do
           printf "%s\\t%s\\n" $f $(samtools idxstats $f | awk -F'\\t' '$1~"^chr[1-9]" {{sum+=$3}} END {{print sum}}') >> {output}
        done
        """

# rule CountJunctionsPerSample:
#     input:

rule QualimapRnaseq:
    input:
        gtf="/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        bam="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        "QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt"
    log:
        "logs/QualimapRnaseq/{sample}.log"
    conda:
        "../envs/qualimap.yml"
    params:
        extra = "-p strand-specific-reverse"
    resources:
        mem_mb = 16000
    shell:
        """
        unset DISPLAY
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} {params.extra} --java-mem-size=12G -outdir QC/QualimapRnaseq/{wildcards.sample}/ &> {log}
        """



use rule featurecounts as featurecounts_chRNA with:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= chRNA_samples.index.unique()),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai", sample= chRNA_samples.index.unique()),
        Comprehensive_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    log:
        "logs/featureCounts_chRNA.log"
    output:
        "featureCounts/Counts.chRNA.txt"
    params:
        extra = "-s 2 -p"

rule MultiQC:
    input:
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= chRNA_samples.index.unique()),
        expand("QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt", sample=chRNA_samples.index),
        expand("FastqFastp/{sample}.fastp.json", sample=chRNA_samples.index)
    log: "logs/Multiqc.log"
    output:
        "Multiqc/multiqc_report.html"
    shell:
        """
        multiqc -f -o Multiqc/ Alignments/STAR_Align/ QC/QualimapRnaseq/ FastqFastp/ &> {log}
        """

rule sampleFragmentLengths:
    input:
        "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "FragLenths/{sample}.1M.test.txt.gz"
    shell:
        """
        samtools view -s 0.1 {input} |  awk '$6 !~ "N" && $9>0 && $9 <1000 {{print $9}}' | shuf -n 1000000 | gzip - > {output}
        """

rule MakeBigwigs_NormalizedToEdgeRFeatureCounts:
    """
    Scale bigwig to base coverage per billion chromosomal reads
    """
    input:
        fai = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
        NormFactorsFile = "../output/QC/ReadCountsAndJunctionsPerSamples.tsv"
    params:
        GenomeCovArgs="-split",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
        BamToBigwigScript = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/scripts/GenometracksByGenotype/BamToBigwig.sh"
        # Region = ""
    # wildcard_constraints:
    #     Phenotype = "|".join(RNASeqPhenotypes)
    shadow: "shallow"
    output:
        bw = "bigwigs/unstranded/{sample}.bw",
        bw_minus = []
    log:
        "logs/MakeBigwigs_unstranded/{sample}.log"
    resources:
        mem = much_more_mem_after_first_attempt
    shell:
        """
        ScaleFactor=$(bc <<< "scale=3;1000000000/$(grep '{input.bam}' {input.NormFactorsFile} | awk 'NR==1 {{print $2}}')")
        {params.BamToBigwigScript} {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs} -scale ${{ScaleFactor}}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """

use rule MakeBigwigs_NormalizedToEdgeRFeatureCounts as MakeBigwigs_NormalizedToEdgeRFeatureCounts_stranded with:
    output:
        bw = "bigwigs/stranded/{sample}.minus.bw",
        bw_minus = "bigwigs/stranded/{sample}.plus.bw"
    log:
        "logs/MakeBigwigs_stranded/{sample}.log"

rule GatherBigwigs:
    input:
        expand("bigwigs/unstranded/{sample}.bw", sample=AllSamples_202310),
        expand("bigwigs/stranded/{sample}.minus.bw", sample=AllNEBNextSamples),
        expand("bigwigs/stranded/{sample}.plus.bw", sample=AllNEBNextSamples)

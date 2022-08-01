
rule GatherChRNAResults:
    input:
        # expand("FastqFastp/{sample}.R1.fastq.gz", sample=chRNA_samples.index),
        # expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=chRNA_samples.index),
        expand("FragLenths/{sample}.1M.test.txt.gz", sample=pd.concat([titration_series_samples, chRNA_samples]).index),
        "Multiqc/multiqc_report.html",
        "featureCounts/Counts.chRNA.txt",
        "../output/QC/ReadCountsAndJunctionsPerSamples.tsv",
        expand("QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt", sample=chRNA_samples.index)

use rule CopyAndMergeFastq as CopyAndMergeFastq_chRNA with:
    input:
        R1 = lambda wildcards: chRNA_samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: chRNA_samples.loc[wildcards.sample]['R2'],
    output:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    wildcard_constraints:
        sample = "|".join(chRNA_samples.index)

rule CountReadsPerSample:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=chRNA_samples.index),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=chRNA_samples.index),
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

# rule BamToBigwig:

# rule bamtobigiwg:
#     input:

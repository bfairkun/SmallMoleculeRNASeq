
rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = lambda wildcards: samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: samples.loc[wildcards.sample]['R2']
    output:
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        R2 = "FastqFastp/{sample}.R2.fastq.gz",
        html = "FastqFastp/{sample}.fastp.html",
        json = "FastqFastp/{sample}.fastp.json"
    params:
    wildcard_constraints:
        sample = "|".join(samples.index)
    resources:
        mem_mb = 8000
    log:
        "logs/fastp/{sample}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2}  -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} &> {log}
        """

rule STAR_Align:
    input:
        index = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        R2 = "FastqFastp/{sample}.R2.fastq.gz"
    output:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        align_log = "Alignments/STAR_Align/{sample}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align_WASP/{sample}.log"
    params:
        GetSTARIndexDir = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/",
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
    resources:
        cpus_per_node = 9,
        mem = 48000,
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix Alignments/STAR_Align/{wildcards.sample}/ --genomeDir {params.GetSTARIndexDir} --readFilesIn {input.R1} {input.R2}  --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 --limitBAMsortRAM 8000000000 {params.ENCODE_params} --outSAMstrandField intronMotif  &> {log}
        """

rule indexBam:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/indexBam/{sample}.log"
    output:
        bai = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    shell: "samtools index {input} &> {log}"

rule ExtractJuncs:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        junc = temp(expand("SplicingAnalysis/leafcutter/juncfiles/chr{chrom}/{{sample}}.junc", chrom=autosomes)),
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",
    params:
        # strand = GetLibStrandForRegtools
        strand = 0
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{sample}.log"
    shell:
        """
        for chrom in {autosomes}
        do
            (regtools junctions extract -m 20 -s {params.strand} -r chr${{chrom}} {input.bam} > SplicingAnalysis/leafcutter/juncfiles/chr${{chrom}}/{wildcards.sample}.junc ) &> {log}
        done
        cat {output.junc} > {output.junc_autosomes}
        """


rule make_leafcutter_juncfile:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=samples.index),
    output:
        "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    params:
        SamplesToRemove = ""
    run:
        import os
        if params.SamplesToRemove:
            SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        else:
            SamplesToRemove=[]
        with open(output[0], "w") as out:
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

use rule make_leafcutter_juncfile as make_leafcutter_juncfile_all with:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=AllSamples),
    output:
        "SplicingAnalysis/leafcutter_all_samples/juncfilelist.txt"

rule leafcutter_cluster:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=samples.index),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_cluster/autosomes.log"
    params:
        "-r SplicingAnalysis/leafcutter/clustering/autosomes/"
    shell:
        """
        python scripts/leafcutter_cluster_regtools_py3.py -j {input.juncfile_list} {params} &> {log}
        """


use rule leafcutter_cluster as leafcutter_cluster_all with:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=AllSamples),
        juncfile_list = "SplicingAnalysis/leafcutter_all_samples/juncfilelist.txt"
    output:
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind_numers.counts.gz"
    params:
        "-r SplicingAnalysis/leafcutter_all_samples/"
    log:
        "logs/leafcutter_cluster_all.log"
    output:
        "SplicingAnalysis/leafcutter_all_samples/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter_all_samples/autosomes/leafcutter_perind_numers.counts.gz"

rule MakeGroupsFiles:
    output:
        "SplicingAnalysis/leafcutter/groupsfiles/{treatment}.txt"
    wildcard_constraints:
        treatment = "|".join(NonControlTreatments)
    run:
        with open(output[0], 'w') as fh:
            TreatmentSamples = samples.loc[samples['Treatment']==wildcards.treatment].index
            for sample in TreatmentSamples:
                _ = fh.write(f"{sample}\t{wildcards.treatment}\n")
            for sample in samples.loc[samples['Treatment']=='DMSO'].index:
                _ = fh.write(f"{sample}\tDMSO\n")


rule annotate_juncfiles:
    input:
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        Comprehensive_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        basic_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=samples.index),
    output:
        junc_file_concat = "SplicingAnalysis/leafcutter/JuncfilesMerged.junc",
        basic = "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.gz",
        comprehensive = "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.comprehensive.bed.gz"
    log:
        "logs/annotate_juncfiles.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        cat {input.juncs} | sort | uniq | bedtools sort -i - > {output.junc_file_concat}
        (regtools junctions annotate {output.junc_file_concat} {input.fa} {input.basic_gtf} | gzip - > {output.basic} ) &> {log}
        (regtools junctions annotate {output.junc_file_concat} {input.fa} {input.Comprehensive_gtf} | gzip - > {output.comprehensive} ) &>> log
        """

rule collapse_juncfiles_allsamples:
    """
    In preparation for regtools annotate using a juncfile with just one line for each unique junction among all junctions discovered across all samples
    """
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=AllSamples),
    output:
        "SplicingAnalysis/FullSpliceSiteAnnotations/ALL_SAMPLES.junc"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Collapse_Juncfiles.R
        """

rule annotate_juncfiles_allsamples:
    input:
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        Comprehensive_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        basic_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        juncs = "SplicingAnalysis/FullSpliceSiteAnnotations/ALL_SAMPLES.junc"
    output:
        basic = "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz",
        comprehensive = "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.comprehensive.bed.gz"
    log:
        "logs/annotate_juncfiles_allsamples.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.juncs} {input.fa} {input.basic_gtf} | gzip - > {output.basic} ) &> {log}
        (regtools junctions annotate {input.juncs} {input.fa} {input.Comprehensive_gtf} | gzip - > {output.comprehensive} ) &>> log
        """

rule Get5ssSeqs:
    """
    Filtered out entries with N in sequence
    """
    input:
        basic = "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.gz",
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
    output:
        temp("SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab")
    shell:
        """
        zcat {input.basic} | awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $1, $2, $3, $1"_"$2"_"$3"_"$6, ".", $6}}' | sort -u | awk -v OFS='\\t' -F'\\t'  '$6=="+" {{$2=$2-4; $3=$2+11; print $0}} $6=="-" {{$3=$3+3; $2=$3-11; print $0}}' | bedtools getfasta -tab -bed - -s -name -fi {input.fa} | grep -v 'N' > {output}
        """

use rule Get5ssSeqs as Get5ssSeqs_allsamples with:
    input:
        basic = "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz",
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
    output:
        temp("SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab")

rule Score5ss:
    input:
        "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab"
    output:
        Tab = "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        WeblogoImg = "../docs/assets/5ssPWM.png"
    conda:
        "../envs/biopython.yml"
    shell:
        """
        python scripts/ScorePWM.py {input} {output.Tab} {output.WeblogoImg}
        """

use rule Score5ss as Score5ss_allsamples with:
    input:
        "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab"
    output:
        Tab = "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        WeblogoImg = "../docs/assets/5ssPWM.allsamples..png"
    conda:
        "../envs/biopython.yml"

rule leafcutter_ds:
    input:
        groupfile = "SplicingAnalysis/leafcutter/groupsfiles/{treatment}.txt",
        numers = "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz",
        Rscript = "scripts/leafcutter/scripts/leafcutter_ds.R"
    output:
        "SplicingAnalysis/leafcutter/differential_splicing/{treatment}_effect_sizes.txt",
        "SplicingAnalysis/leafcutter/differential_splicing/{treatment}_cluster_significance.txt"
    threads: 4
    resources:
        ntasks = 5
    params:
        "-i 3 -g 3"
    log:
        "logs/leafcutter_ds/{treatment}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript {input.Rscript} -p {threads} -o  SplicingAnalysis/leafcutter/differential_splicing/{wildcards.treatment} {params} {input.numers} {input.groupfile} &> {log}
        """

rule featurecounts:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= samples.index),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai", sample= samples.index),
        Comprehensive_gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "featureCounts/Counts.txt"
    threads:
        8
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts.log"
    params:
        extra = ""
    shell:
        """
        featureCounts {params.extra} -T {threads} --ignoreDup --primary -a {input.Comprehensive_gtf} -o {output} {input.bam} &> {log}
        """

rule DE_testing:
    input:
        "featureCounts/Counts.txt"
    output:
        results = "DE_testing/Results.txt.gz",
        counts = "DE_testing/Counts.mat.txt.gz"
    log:
        "logs/DE_testing.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/DE_edgeR.R {input} {output.results} {output.counts} &> {log}
        """

rule Bam_list:
    input:
        bam = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= samples.index),
        bai = expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai", sample= samples.index),
    output:
        "../data/bam_list.txt"
    shell:
        """
        ls -1 Alignments/STAR_Align/*/Aligned.sortedByCoord.out.bam | awk -v OFS='\\t' '{split($1, a, "/"); print $1, a[3]}' > {output}
        """

rule JuncFilesToBed:
    input:
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",
    output:
        bed = "SplicingAnalysis/juncbedfiles/{sample}.bed.gz"
    shell:
        """
        cat <(echo 'track name="{wildcards.sample}"  itemRgb="On" graphType=junctions') {input} | gzip - > {output}
        """

rule MorePermissiveLeafcutterClustering:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{sample}.junc",sample=samples.index),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering_permissive/autosomes/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering_permissive/autosomes/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_cluster/autosomes.log"
    shell:
        """
        python scripts/leafcutter_cluster_regtools_py3.py -m 5 -j {input.juncfile_list} -r SplicingAnalysis/leafcutter/clustering_permissive/autosomes/ &> {log}
        """

rule spliceq:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        gtf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf.ensembl",
    output:
        "SplicingAnalysis/SpliceQ/ByChrom/{chrom}/{sample}.txt"
    conda:
        "../envs/spliceq.yaml"
    log:
        "logs/spliceq/{sample}.{chrom}.log"
    shell:
        """
        SPLICE-q.py -b {input.bam} -g {input.gtf} -o {output} -p 1  -c 0 -i -x {wildcards.chrom}
        """

rule GatherSpliceq:
    input:
        expand("SplicingAnalysis/SpliceQ/ByChrom/{chrom}/{{sample}}.txt", chrom=range(1,23))
    output:
        "SplicingAnalysis/SpliceQ/Merged/{sample}.txt.gz"
    shell:
        """
        awk 'NR==1 || FNR>1' {input} | gzip - > {output}
        """

rule MergeSpliceQ_Quantifications:
    input:
        expand("SplicingAnalysis/SpliceQ/Merged/{sample}.txt.gz", sample=AllSamples),
    output:
        "SplicingAnalysis/SpliceQ/MergedTable.txt.gz"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rsript 
        """

rule TidyDoseResponseData:
    input:
        "bigwigs/BigwigList.tsv",
        "bigwigs/BigwigList.groups.tsv",
        "SplicingAnalysis/leafcutter_all_samples/PSI.table.bed.gz",
        "featureCounts/Counts.titration_series.txt",
        "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        "SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz",
    output:
        "DoseResponseData/LCL/TidySplicingDoseData.txt.gz",
        "DoseResponseData/LCL/TidyExpressionDoseData.txt.gz",
    log:
        "logs/TidyDoseResponseData.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/TidyAndSpearmanCorTitrationData.R &> {log}
        """

## Check how ClinVar splice site annotations relate to splice sites:
# zcat ClinVar/PangolinResults.tsv.gz | grep 'splice_donor' | awk -v OFS='\t' '{print "chr"$1, $2, $2+1}' | sort -u | bedtools sort -i - | bedtools closest -a - -b <( bedtools flank -g /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai  -i <(zcat SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.gz | awk -v OFS='\t' 'NR>1 {print $1,$2,$3,".",".",$6}' ) -l 1 -r 0 | sort -u | bedtools sort -i -  ) -D b | less -S

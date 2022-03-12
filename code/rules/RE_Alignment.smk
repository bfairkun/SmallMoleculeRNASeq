
def get_FilterBam_Inputs(wildcards):
    sample_name = samples.query('sample == @wildcards.sample').index[0]
    return expand("Alignments/STAR_Align/{S}/Aligned.sortedByCoord.out.bam", S=sample_name)

rule FilterBam_NoWASP:
    message: '''### Filter BAM on flag, mapq, chrom ### '''
    input: get_FilterBam_Inputs
    output: 
        bam = temp("RE_Results/Alignments/filter/{sample}.bam")
    log: "logs/FilterBam_NoWASP/{sample}.log"
    params:
        inbam = "RE_Results/Alignments/{sample}.bam",
        chroms = CHROMS
    threads: 4
    resources: time=500, mem_mb=12000, cpu=4
    shell:
        '''
        samtools view -b -h -@ {threads} \
            -F 256 -F 2048 -q 5 \
            -o {output.bam}.temp.bam \
            {input} &> {log}
        samtools index -@ {threads} {output.bam}.temp.bam &>> {log}        
        samtools view -b -h -@ {threads} -o {output} {output}.temp.bam {params.chroms} &>> {log}
        samtools index -@ {threads} {output.bam}  &>> {log}
        rm {output.bam}.temp* &>> {log}
        '''



# This is needed for alignments done by Ben
# Not needed for mRNA alignments done by me
rule AddReadGroupTags:
    message: '### Add RG tags'
    input:
        bam = "RE_Results/Alignments/filter/{sample}.bam"
    output: 
        bam = temp("RE_Results/Alignments/AddRG/{sample}.bam"),
        bai = temp("RE_Results/Alignments/AddRG/{sample}.bam.bai")
    log: "logs/AddReadGroupTags/{sample}.log"
    params:
        RGID = '{sample}',
        RGSM = '{sample}',
        RGPL = "ILLUMINA",
        RGLB = '{sample}',
        RGPU = '{sample}' 
    threads: 1
    resources: time = 500, mem_mb = 15000, cpu = 1
    shell:
        '''
        gatk AddOrReplaceReadGroups -I {input.bam} -O {output.bam} \
            --RGID {params.RGID} --RGSM {params.RGSM} --RGPL {params.RGPL} --RGLB {params.RGLB} --RGPU {params.RGPU}
        samtools index {output}
        '''
    

# mark duplicates
rule MarkDups:
    message: '''### MarkDuplications ###'''
    input:
        bam = "RE_Results/Alignments/AddRG/{sample}.bam",
    output:
        bam = temp("RE_Results/Alignments/MarkDups/{sample}_markDups.bam"),
        metrics = "RE_Results/Alignments/MarkDups/{sample}_MarkDupsMetrics.txt"
    log: "logs/MarkDups/{sample}.log"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    threads: 4
    resources: time=2000, mem_mb=25000, cpu=4
    shell:
        '''
        gatk MarkDuplicatesSpark  \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --tmp-dir {params.tmp} \
            --remove-sequencing-duplicates &> {log}
        '''

# # BQSR on RNA-seq requires first split reads that are mapped to multiple loci, 
# # indicated by the N number of cigar string. 
# # SplitNCigarReads
rule SplitNCigarReads:
    message: '''### Split chimeric RNA-seq reads using the N Cigar string ###'''
    input:
        bam = rules.MarkDups.output.bam,
        ref = config['FA_HS38']
    output:
        bam = temp("RE_Results/Alignments/SplitNCigar/{sample}.bam")
    log: "logs/SplitNCigarReads/{sample}.log"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    threads: 1
    resources: time=2000, mem_mb=25000, cpu=1
    shell:
        '''
        gatk SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam} \
            --tmp-dir {params.tmp} &> {log}
        '''

# BQSR step 1
# Note GATK BaseRecalibrator will not take in any reads that has N in cigar
rule BaseRecalibrator:
    message: '''### Compute covariate matrix of Base Recalibration ###'''
    input:
        bam = rules.SplitNCigarReads.output.bam,
        ref = config['FA_HS38'],
        known_sites = config['dbSNP']
    output:
        recal_file = "RE_Results/BQSR/{sample}_covariates.tab"
    log: "logs/BaseRecalibrator/{sample}.log"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory errorthreads: 1
    resources: time = 2000, mem_mb = 25000, cpu=1
    shell:       
        '''
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.known_sites} \
            --tmp-dir {params.tmp} \
            -O {output.recal_file} &> {log}
        '''

# BQSR step 2
# do not use temp() function on output. Remove bams manually after merged.
rule ApplyBQSR:
    message: '''### Apply Base Recalibration, output calibrated BAM ###'''
    input:
        bam = rules.SplitNCigarReads.output.bam,
        ref = config['FA_HS38'],
        recal_file = rules.BaseRecalibrator.output.recal_file
    output:
        bam = "RE_Results/BQSR/{sample}_recal.bam"
    log: "logs/ApplyBQSR/{sample}.log"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory errorthreads: 1
    resources: time=2000, mem_mb=25000, cpu=1
    shell:
        '''
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.recal_file} \
            --tmp-dir {params.tmp} \
            -O {output.bam} &> {log}
        '''

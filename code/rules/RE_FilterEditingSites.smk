

rule Gather1stPassSites:
    """
    Gather all the editing sites from first pass.
    This results in a VCF file that's compatable for STAR wasp.
    Only require ALT depth > 0 and AF > 0; 
    """
    message: '''### Make VCF file from A>G editing sites'''
    input: expand("RE_Results/BCFCall/FinalAnno/{sample}.bed", sample=SAMPLE_NAMES)
    output: "RE_Results/Remap/STAR-remap-variant.vcf"
    log: 'logs/Gather1stPassSites.log'
    threads: 1
    resources: time=120, mem_mb=8000, cpu=1
    script: "../scripts/makeA2G_VCF.R"



# This step requires a lot of tmp space when bam files are large, e.g. with chRNA
rule BAM_to_Fastq:
    message: '''### Extract paired end fastq from BAM'''
    input: 'RE_Results/BQSR/{sample}_recal.bam'
    output: 
        R1 = temp("RE_Results/BamToFastq/{sample}.R1.b2f.fastq.gz"),
        R2 = temp("RE_Results/BamToFastq/{sample}.R2.b2f.fastq.gz"),
    log: 'logs/BAM_to_Fastq/{sample}.log'
    params:
        collate_tmp = "/scratch/midway2/chaodai/TMP/{sample}.collate.tmp.bam",
        prefix = "/scratch/midway2/chaodai/TMP/{sample}"
    threads: 4
    resources: time=2000, mem_mb=20000, cpu=4
    shell:
        '''
        samtools collate -@ {threads} -o {params.collate_tmp} {input} {params.prefix} &> {log} 
        samtools fastq -@ {threads} -1 {output.R1} -2 {output.R2} -0 /dev/null -s /dev/null {params.collate_tmp} &>> {log}
        rm {params.collate_tmp}
        '''


rule Remap_STAR:
    """
    Remap reads that overlap 1st pass A>G sites
    """
    message: '''### Remap using STARS ###'''
    input:
        R1 = "RE_Results/BamToFastq/{sample}.R1.b2f.fastq.gz",
        R2 = "RE_Results/BamToFastq/{sample}.R2.b2f.fastq.gz",
        A2G_VCF = "RE_Results/Remap/STAR-remap-variant.vcf"
    output: temp("RE_Results/Remap/{sample}.bam")
    log: 'logs/Remap_STAR/{sample}.log'
    params:
        GENOMEDIR = config['STAR_INDEX_DIR'],
        REF_FA = config['FA_HS38'],
        OUT_PREFIX = "RE_Results/Remap/{sample}.",
        read_files_comm = "zcat",
        limitOutSJcollapsed = 5000000, # adding this for nuclear RNA, otherwise report error sometimes
    threads: 12
    resources: time=2100, mem_mb=40000, cpu=12, partition="bigmem2"
    shell:
        '''
        STAR="/software/STAR-2.7.7a-el7-x86_64/bin/STAR"
        $STAR \
            --runThreadN {threads} \
            --genomeDir {params.GENOMEDIR} \
            --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand {params.read_files_comm} \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD RG vW \
            --outSAMunmapped None \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.OUT_PREFIX} \
            --limitOutSJcollapsed {params.limitOutSJcollapsed} \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} PL:ILLUMINA LB:{wildcards.sample} PU:{wildcards.sample} \
            --varVCFfile {input.A2G_VCF} \
            --waspOutputMode SAMtag &> {log}
        
        sleep 10
        samtools sort -@ {threads} -o {output} {params.OUT_PREFIX}Aligned.out.bam &>> {log}
        sleep 60
        samtools index -@ {threads} {output} &>> {log}
        sleep 10
        rm {params.OUT_PREFIX}Aligned.out.bam &>> {log}
        
        '''


### filter remapped reads to include only reads that pass WASP
rule Filter_Remapped_Bam:
    message: '### Remove alignments failing to pass WASP tag in the 2nd STAR mapping'
    input: "RE_Results/Remap/{sample}.bam"
    output: "RE_Results/Remap/filterBam/{sample}.bam"
    log: 'logs/Filter_Remapped_Bam/{sample}.log'
    threads: 4
    resources: time=60, mem_mb=20000, cpu=4, partition="broadwl"
    shell:
        '''
        samtools view -@ {threads} -h -F 256 -F 2048 -q 20 \
          -e " [vW] == 1" \
          -o {output} {input} &> {log}
        sleep 10
        samtools index -@ {threads} {output} &>> {log}
        '''

### filter editing sites based on remap consistency 
rule Filter_Remap_Consistency:
    message: '### Filter sites based on remap requirements'
    input:
        BED = "RE_Results/BCFCall/FinalAnno/{sample}.bed",
        BAMPRE = 'RE_Results/BQSR/{sample}_recal.bam',
        BAMPOST = "RE_Results/Remap/filterBam/{sample}.bam"
    output:
        tab = "RE_Results/RemapConsistency/{sample}.txt",
        bed = "RE_Results/RemapConsistency/{sample}.bed"
    log: 'logs/Remap_Consistency/{sample}.log'
    threads: 1
    resources: time=100, mem_mb=30000, cpu=1, partition="broadwl"
    shell:
        '''
        python scripts/checkRemap.py \
            --BED {input.BED} --bamPre {input.BAMPRE} --bamPost {input.BAMPOST} \
            --outTab {output.tab} --outBED {output.bed} &> {log}
        '''

rule Filter_Primer_Error:
    message: '### Remove sites supported primer errors (6bp with 5 prim)'
    input: 
        BED = "RE_Results/RemapConsistency/{sample}.bed",
        BAM = "RE_Results/Remap/filterBam/{sample}.bam"
    output: 'RE_Results/removePrimerError/{sample}.bed'
    log: 'logs/Primer_Error/{sample}.log'
    threads: 1
    resources: time=100, mem_mb=15000
    shell:
        '''
        python scripts/check5Preads.py \
            --BED {input.BED} --bam {input.BAM} --th 0.2 --outBED {output} &> {log}
        '''

rule Filter_Simple_Repeats:
    message: '### Remove simple repeats'
    input: 'RE_Results/removePrimerError/{sample}.bed'
    output: 'RE_Results/removeSimpleRepeats/{sample}.bed' 
    log: 'logs/Simple_Repeats/{sample}.log'
    params: 
        SimpleRepeats = config['SIMPLE_REPEATS']
    threads: 1
    shell:
        '''
        intersectBed -wa -v -a {input} -b {params.SimpleRepeats} 1> {output} 2>{log}
        '''

# merge SJ tab files and produce bed for later use.
rule GatherSpliceJunctions:
    message: '''### Make Bed files of coordinates of 4bp adjacent to splice sites in introns '''
    input: "Alignments/STAR_Align/{sample}/SJ.out.tab"
    output: "RE_Results/Alignment/STAR/SJ/{sample}.SJ.4bpIntron.bed"
    log: 'logs/GatherSpliceJunctions/{sample}.log'
    threads: 1
    shell:
        '''
        Rscript scripts/SJ_v2.R {input} {output} &> {log}
        '''


rule Filter_Intron_SJ:
    """
    Filter intronic sites within 4bp of splice sites, 
    using specific intron splice sites for each sample produced by STAR
    """
    message: '### Remove intronic sites with 4bp of splice sites'
    input: 
        Bed = 'RE_Results/removeSimpleRepeats/{sample}.bed',
        IntronSJ = "RE_Results/Alignment/STAR/SJ/{sample}.SJ.4bpIntron.bed"
    output: 'RE_Results/removeIntronSJ/{sample}.bed'
    log: 'logs/Filter_Intron_SJ/{sample}.log'
    threads: 1
    shell:
        '''
        intersectBed -wa -v -a {input.Bed} -b {input.IntronSJ} 1> {output} 2>{log}
        '''    

### filter homopolymers >= 5
rule Filter_homopolymers:
    message: '### Remove homopolymers of length >=5 '
    input: 'RE_Results/removeIntronSJ/{sample}.bed'
    output: 'RE_Results/removeHomopolymers/{sample}.bed'
    log: 'logs/Filter_Homopolymers/{sample}.bed'
    params: 
        Homopolymer = config['HOMOPOLYMER']
    threads: 1
    shell:
        '''
        intersectBed -wa -v -a {input} -b {params.Homopolymer} > {output}
        '''


rule bcf_to_bed:
    """
    Convert bcf to bed for selected directories. Note converting VCF, 1-based, into BED, 0-based format
    """
    message: '### Convert bcf to bed for some directories'
    input: "RE_Results/BCFCall/raw/{sample}.bcf"
    output: 'RE_Results/BCFCall/raw/{sample}.bed'
    log: 'logs/bcf_to_bed/{sample}.log'
    threads: 1
    shell:
        '''
            bcftools view -H --threads {threads} {input} 2>{log} | \
                awk '   BEGIN {{ OFS="\t" }};
                        {{  
                            start=$2-1;
                            end=$2;
                            print $1, start, end, $4 "," $5, int($6), "."
                        }}
                    ' 1> {output} 2>{log}
        '''

rule CountCoverage:
    message: '### Count coverage at editing site'
    input: 
        bed = 'RE_Results/removeHomopolymers/{sample}.bed',
        bam = 'RE_Results/Remap/filterBam/{sample}.bam'
    output: 'RE_Results/countCoverage/{sample}.txt'
    log: 'logs/CountCoverage/{sample}.log'
    threads: 1
    shell:
        '''
            python scripts/countCoverage.py \
                --BAM {input.bam} --BED {input.bed} --outTab {output} &>{log}
        '''

rule QuantEditingLevel:
    message: '### Quantify editing level'
    input: 
        bed = 'RE_Results/BCFCall/raw/{sample}.bed',
        cnt = 'RE_Results/countCoverage/{sample}.txt'
    output: 'RE_Results/EditLevel/{sample}.txt'
    log: 'logs/QuantEditingLevel/{sample}.log'
    script: '../scripts/quantEditLevel.R'
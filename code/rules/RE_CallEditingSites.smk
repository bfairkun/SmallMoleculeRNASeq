# call variants
rule CallVariant:
    message: '''### Call variants ###'''
    input:
        bam = 'RE_Results/BQSR/{sample}_recal.bam',
        ref = config['FA_HS38']
    output: 
        bcf = "RE_Results/BCFCall/raw/{sample}.bcf"
    log: 'logs/CallVariant/{sample}.log'
    params: region = "" # use "-r chr21" for limited regions
    threads: 6
    resources: time=2000, mem_mb=20000, cpu=6
    shell:
        '''
        bcftools mpileup -q 20 -Q 20 {params.region} --threads {threads} \
            -a FORMAT/AD,FORMAT/SP,INFO/AD \
            -Ou -f {input.ref} {input.bam} 2> {log} | \
            bcftools call -m -v -Ob -o {output.bcf} --threads {threads} &>> {log}
        bcftools index --threads {threads} {output.bcf} &>> {log}
        '''


CHROMS_str = ",".join(CHROMS)
rule Filter1:
    """
    filter variants using these criteria: 
        - 1. SNPs only (remove indels)
        - 2. A to G or T to C variants
        - DP (total read depth) > 9, 
        - within chr1 - chrY
        - QUAL >= 20
    """
    message: '''### Apply 1st set of filters ###'''
    input:
        bcf = "RE_Results/BCFCall/raw/{sample}.bcf"
    output:
        "RE_Results/BCFCall/filter1/{sample}.bcf"
    log: "logs/Filter1/{sample}.log"
    params:
        chrom = CHROMS_str
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    shell:
        '''
        bcftools filter --threads {threads} -r {params.chrom} \
            -i 'INDEL=0 && QUAL > 19 && DP>9 && ((REF="A" && ALT="G") || (REF="T" && ALT="C"))' \
            -Ob -o {output} {input.bcf} &> {log}
        bcftools index --threads {threads} {output} &>> {log}
        '''


rule Ann_1KGP3:
    """
    Annotate with dbSNP and remove any COMMON=1 variants
    ID=COMMON,Number=1,Type=Integer,Description="RS is a common SNP.  
    A common SNP is one that has at least one 1000Genomes population with a minor 
    allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.  
    """
    message: '''### Annotate using 1KGP vcf, remove common variants (AF >1%) ###'''
    input: 
        bcf = "RE_Results/BCFCall/filter1/{sample}.bcf",
        ann = config['KGP3']
    output: "RE_Results/BCFCall/Ann1KGP/{sample}.bcf"
    log: "logs/Ann_1KGP3/{sample}.log"
    params:
        temp = "RE_Results/BCFCall/Ann1KGP/{sample}.temp.bcf"
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    shell:
        '''
        bcftools annotate -a {input.ann} -c ID,+INFO {input.bcf} 2> {log} | \
            bcftools filter -i 'ID="."' -Ob -o {output} &>> {log}
        bcftools index {output} &>> {log}
        '''


rule CheckDB:
    """
    Check existing RNA-editing database. Here I used REDIportal database
    Columns after bcftools query: chr, pos, id, ref, alt, qual, dp, ad, dbname, mq
    Columns after awk: chr, bedStart, bedEnd, id, score, strand, ref, alt, qual, dp, AD of ref, AD of alt, AF, dbname, mq
    """
    message: '''### Check existing RNA Editing database ###'''
    input: 
        bcf = "RE_Results/BCFCall/Ann1KGP/{sample}.bcf",
        db = config['DB'],
        header = config['DB_HEADER']
    output: 
        bcf = "RE_Results/BCFCall/checkDB/{sample}.bcf",
        tab = "RE_Results/BCFCall/checkDB/{sample}.bed" #0 indexed
    log: "logs/CheckDB/{sample}.log"
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    shell:
        '''
        bcftools annotate -a {input.db} -h {input.header} \
            -c CHROM,FROM,TO,INFO/DBNAME,INFO/BEDSCORE,INFO/STRAND \
            -Ob -o {output.bcf} {input.bcf} &> {log}
        bcftools index {output.bcf} &>> {log}
        
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%DP\t%AD\t%DBNAME\t%MQ\n' {output.bcf} 2>> {log} | \
            awk 'BEGIN {{OFS="\t"}}; {{
                        start=$2-1; 
                        end=$2; 
                        split($8, AD, ","); 
                        AF=AD[2]/(AD[1]+AD[2]); 
                        print $1,start, end, $3,999, ".",$4,$5,$6,$7,AD[1],AD[2],AF,$9,$10
                        }} ' 1> {output.tab} 2>> {log}
       '''

rule Filter2:
    """
    Additional filters. e.g. mininum AD >= 1
    outfile column names: "chr", "start", "end", "id", "score", "strand", 
            "ref", "alt", "qual", "dp", "ad_ref", "ad_alt", "af", "rnadb", "mapq"
    """
    message: '''### Apply 2nd set of filters ###'''
    input: 'RE_Results/BCFCall/checkDB/{sample}.bed' # bed file, so 0 index
    output: "RE_Results/BCFCall/Filter2/{sample}.bed" # bed, 0 index
    log: 'logs/Filter2/{sample}.log'
    threads: 1
    script:
        "../scripts/FinalRNASites_bcf.R"


rule Gene_Features:
    """
    Columns before 1st intersectBed: "chr", "start", "end", "id", "score", "strand", "ref", "alt"
    Columns after cut -f 1-8,12,15 "chr", "start", "end", "id", "score", "strand", "ref", "alt", "geneName", "featureName"
    Columns after 2nd intersectBed: "chr", "start", "end", "id", "score", "strand", "ref", "alt", "geneName", "featureName", "repeatName"
    Columns after params.R: "chr", "start", "end", "id", "score", "strand", "ref", "alt", "geneName", "featureName", "repeatName", "repeatFamily"
    """
    message: '''### Add gene features and output tab file ###'''
    input:
        bed = 'RE_Results/BCFCall/Filter2/{sample}.bed', # bed file 0 indexed
        gencode = config['GENCODE'],
        repeats = config['REPEATS'],
        repeatsLookup = config['REPEATS_LOOKUP']
    output: "RE_Results/BCFCall/FinalAnno/{sample}.bed"
    log: "logs/Gene_Features/{sample}.log"
    params:
        R = "scripts/GeneAnnotation.R",
        igvbed = "RE_Results/BCFCall/FinalAnno/{sample}.igv.bed"
    shell:
        '''
        intersectBed -wao -a {input.bed} -b {input.gencode} 2> {log} | \
            cut -f 1-8,12,15 | uniq | \
            intersectBed -wao -a - -b {input.repeats} 2>> {log} | \
            cut -f 1-10,14 | uniq > {output}
        
        Rscript {params.R} {output} {input.repeatsLookup} &>> {log}
        
        sleep 10
        
        awk 'BEGIN {{OFS="\t"}}; {{name=$4"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12; print $1,$2,$3,name,$5,$6}}' {output} 1> {params.igvbed} 2>> {log}
        '''




rule GatherChemCLIPResults:
    input:
        expand("ChemCLIP/bigwigs/{sample}.bw", sample=["Treatment", "Control"])


rule DownloadChemCLIPFastq:
    input:
    output:
        Treatment = "ChemCLIP/Fastq/Treatment.fastq.gz",
        Control = "ChemCLIP/Fastq/Control.fastq.gz"
    log:
        "logs/DownloadChemCLIPFastq.log"
    shell:
        """
        wget -O- "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/009/SRR6364169/SRR6364169.fastq.gz" > {output.Control}
        wget -O- "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR636/000/SRR6364170/SRR6364170.fastq.gz" > {output.Treatment}
        """

rule STAR_Align_ChemCLIP:
    input:
        index = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "ChemCLIP/Fastq/{sample}.fastq.gz",
    output:
        bam = "ChemCLIP/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        align_log = "ChemCLIP/STAR_Align/{sample}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align_ChemCLIP/{sample}.log"
    params:
        GetSTARIndexDir = "/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex",
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
    resources:
        cpus_per_node = 9,
        mem = 48000,
    wildcard_constraints:
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix ChemCLIP/STAR_Align/{wildcards.sample}/ --genomeDir {params.GetSTARIndexDir} --readFilesIn {input.R1}   --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 --limitBAMsortRAM 8000000000 {params.ENCODE_params} --outSAMstrandField intronMotif  &> {log}
        """

rule indexChemCLIP_bams:
    input:
        bam = "ChemCLIP/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bai = "ChemCLIP/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """

rule ChemCLIP_bigwigs:
    input:
        fai = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = "ChemCLIP/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "ChemCLIP/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        bw = "ChemCLIP/bigwigs/{sample}.bw"
    shell:
        """
        /project2/yangili1/bjf79/GenometracksByGenotype/BamToBigwig.sh {input.fai} {input.bam} {output.bw} -split
        """

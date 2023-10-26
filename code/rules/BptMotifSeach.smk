
rule GetMemeInput_BPRegion:
    input:
        juncs = "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.gz",
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        "Meme/SeqsTab/memetest.50.fa.tab"
    shell:
        """
        zcat {input.juncs} | awk -F'\\t' -v OFS='\\t' 'NR>1 && $14==1 {{$4="."; $5="."; print $0}}' | head -10000 | awk -F'\\t' -v OFS='\\t' '$6=="+" {{$2=$3-51; $3=$3-1; print $0}} $6=="-" {{$3=$2+50; print $0}}' | awk -v OFS='\\t' '{{print $1,$2,$3, "BP", ".", $6}}' |sort | uniq |  bedtools getfasta -bed - -fi {input.fa}  -s -tab > {output}
        """

rule MakeMemeInput:
    input:
        "Meme/SeqsTab/memetest.50.fa.tab"
    output:
        fa = "Meme/Fasta/BP.fa",
        psp = "Meme/PSP/BP.psp"
    conda:
        "../envs/r_essentials.yml"
    params:
        width = 5
    shell:
        """
        Rscript scripts/GetBranchWindow.R {input} {output.psp} {output.fa} {params.width}
        """

rule Meme:
    input:
        fa = "Meme/Fasta/BP.fa",
        psp = "Meme/PSP/BP.psp"
    params:
        width = 5
    output:
        "Meme/results/BP/meme.html"
    conda:
        "../envs/meme.yaml"
    shell:
        """
        meme {input.fa} -w {params.width} -psp {input.psp} -oc Meme/results/BP -dna -nmotifs 1 -mod oops
        """

rule MakeSitesFilesFor5ss:
    input:
        "../output/EC50Estimtes.FromPSI.txt.gz"
    output:
        B = "DonorMotifSearches/Sites/BranaplamSpecific.txt",
        R = "DonorMotifSearches/Sites/RisdiplamSpecific.txt"
    log:
        "logs/MakeSitesFilesFor5ss.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/MakeDonorSitesForMeme.R {input} DonorMotifSearches/Sites/ &> {log}
        """

rule Sites2memeForDonorSites:
    input:
        B = "DonorMotifSearches/Sites/BranaplamSpecific.txt",
        R = "DonorMotifSearches/Sites/RisdiplamSpecific.txt"
    output:
        "DonorMotifSearches/Motifs.meme"
    conda:
        "../envs/meme.yaml"
    shell:
        """
        sites2meme DonorMotifSearches/Sites > {output}
        """

rule fimo_DonorSites:
    input:
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        motifs = "DonorMotifSearches/Motifs.meme"
    output:
        bed = "DonorMotifSearches/DonorHits.fimo.bed.gz",
        tbi = "DonorMotifSearches/DonorHits.fimo.bed.gz.tbi"
    log:
        "logs/fimo_DonorSites.log"
    conda:
        "../envs/meme.yaml"
    shell:
        """
        (fimo --text {input.motifs} {input.fa} | awk -v OFS='\\t' 'NR>1 {{print $2, $3, $4, $1"_"$8, $6, $5}}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed}) &>> {log}
        (tabix -p bed {output.bed}) &>> {log}
        """

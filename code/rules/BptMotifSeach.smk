
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

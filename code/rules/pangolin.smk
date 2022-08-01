rule DownloadVcf:
    output:
        vcf = "ClinVar/clinvar.vcf.gz",
        tbi = "ClinVar/clinvar.vcf.gz.tbi"
    shell:
        """
        wget -O {output.vcf} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
        wget -O {output.tbi} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
        """

rule SplitClinVarVcfByChrom:
    input:
        vcf = "ClinVar/clinvar.vcf.gz",
        tbi = "ClinVar/clinvar.vcf.gz.tbi"
    output:
        vcf = "ClinVar/ByChrom/{chrom}.vcf",
    params:
    resources:
        mem_mb = 8000
    log:
        "logs/SplitClinVarVcfByChrom/{chrom}.log"
    shell:
        """
        (bcftools view {params} -O v {input.vcf} {wildcards.chrom} > {output.vcf} ) &> {log}
        """

rule pangolin:
    input:
        vcf = "ClinVar/ByChrom/{chrom}.vcf",
        fa = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        db = "/project2/yangili1/bjf79/20210618_Pangolin/gencode.v38.basic.annotation.db"
    log:
        "logs/pangolin/{chrom}.log"
    output:
        vcf = "ClinVar/ByChrom/{chrom}.pangolin.vcf"
    shell:
        """
        bash -c '. $HOME/.bashrc; conda activate pangolin'
        sbatch --wait --err {log} --partition=gpu2 --gres=gpu:2 --wrap="set -xe; pangolin {input.vcf} {input.fa} {input.db} ClinVar/ByChrom/{wildcards.chrom}.pangolin"
        """

rule compress_and_index_pangolin_results:
    input:
        original = "ClinVar/ByChrom/{chrom}.vcf",
        pangolin = "ClinVar/ByChrom/{chrom}.pangolin.vcf"
    output:
        vcf = "ClinVar/ByChrom/{chrom}.pangolin.vcf.gz",
        tbi = "ClinVar/ByChrom/{chrom}.pangolin.vcf.gz.tbi",
    shell:
        """
        bcftools view -O z {input.pangolin} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule MergePangolingResults:
    input:
        vcf = expand("ClinVar/ByChrom/{chrom}.pangolin.vcf.gz", chrom=autosomes),
        tbi = expand("ClinVar/ByChrom/{chrom}.pangolin.vcf.gz.tbi", chrom=autosomes),
    output:
        vcf = "ClinVar/clinvar.pangolin.vcf",
    shell:
        """
        bcftools concat {input.vcf} > {output.vcf}
        """


rule ParsePangolinVcf:
    input:
        vcf = "ClinVar/clinvar.pangolin.vcf",
    output:
        "../output/ClinVar.PangolinResults.tsv.gz"
    shell:
        """
        python scripts/ExtractFromVcf.py {input.vcf} {output}
        """


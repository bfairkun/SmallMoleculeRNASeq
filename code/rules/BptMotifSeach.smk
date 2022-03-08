rule GetBpRegions:
    """
    Based on persuing bps originally from Pineda et al, I think a reasonable region to assume bps are located are between 16-42 bp upstream of 3'ss. In particular, filtering for the most believable branches from my prior filter of 5-50 bp, 90% of the branch supporting reads map the branchpoint to this 16-42 region upstream of 3'ss.
    """
    input:
    shell:
        """
        zcat SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.gz | awk -F'\t' -v OFS='\t' '$6=="+" {$2=$3-42; $3=$3-16; print $0} $6=="-" {$3=$2+42; $2=$2+16; print $0}' | bedtools getfasta -bed - -fi /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa -s -tab |  grep 'CTAAC' | sort | uniq | head
        """


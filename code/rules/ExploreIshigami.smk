
rule DownloadIshigamiGEO:
    output:
        A = "Ishigami/GSE221868_SE.MATS.JC.txt.gz",
        B = "Ishigami/GSE221868_elp1_psi.csv.gz",
        C = "Ishigami/GSE221868_exon_skipping.filtered.psi.csv.gz",
        D = "Ishigami/GSE221868_smn2_psi.csv.gz",

    shell:
        """
        wget -O {output.A} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221868/suppl/GSE221868%5FSE%2EMATS%2EJC%2Etxt%2Egz
        wget -O {output.B} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221868/suppl/GSE221868%5Felp1%5Fpsi%2Ecsv%2Egz
        wget -O {output.C} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221868/suppl/GSE221868%5Fexon%5Fskipping%2Efiltered%2Epsi%2Ecsv%2Egz
        wget -O {output.D} https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221868/suppl/GSE221868%5Fsmn2%5Fpsi%2Ecsv%2Egz
        """



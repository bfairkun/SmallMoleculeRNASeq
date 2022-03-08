#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : DE_test
# @created     : Tuesday Jan 04, 2022 10:02:01 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/Counts.txt DE_tests.txt.gz DE_tests.mat.counts.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(edgeR)

f_in <- args[1]
f_out <- args[2]
count_out <- args[3]

rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/(.+?)/Aligned\\.sortedByCoord\\.out\\.bam", "\\1")
    )
}

genes <- read_tsv("../data/Genes.list.txt")

counts <- read_tsv(f_in, comment="#") %>%
    rename_at(vars(-(1:6)),  rename_STAR_alignment_samples) %>%
    select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
    mutate(ensembl_gene_id = str_replace(Geneid, "(.+?)\\.\\d+$", "\\1")) %>%
    left_join(genes, by="ensembl_gene_id") %>%
    mutate(Geneid = paste(Geneid, hgnc_symbol, sep="_")) %>%
    column_to_rownames("Geneid") %>%
    select(-ensembl_gene_id, -hgnc_symbol)

counts %>%
    rownames_to_column("Geneid") %>%
    write_tsv(count_out)

group <- counts %>% colnames() %>%
    str_replace("-\\d+", "") %>% as.factor() %>%
    relevel("DMSO")

design <- model.matrix(~group)

y <- DGEList(counts, group = group)

mean.cpm <- y %>%
    cpm(log=T) %>%
    apply(1, mean)


# filtered <- y[mean.cpm>1, , keep.lib.sizes=FALSE] %>%
#     calcNormFactors()

keep <- filterByExpr(y)
filtered <- y[keep, , keep.lib.sizes=FALSE]


# CONTRASTS <- makeContrasts( C2 = groupC2-groupDMSO,
#                             CUOME = groupCUOME-groupDMSO,
#                             SM2 = groupSM2-groupDMSO,
#                             A10 =groupA10-groupDMSO,
#                             levels = design )

CONTRASTS <- makeContrasts( C2 = groupC2,
                            CUOME = groupCUOME,
                            SM2 = groupSM2,
                            A10 = groupA10,
                            levels = design )

fit <- filtered %>%
    estimateDisp(design) %>%
    glmQLFit(design)

N <- ncol(CONTRASTS)
results <- vector("list", N) %>%
    setNames(colnames(CONTRASTS))

for (i in 1:ncol(CONTRASTS)){
    results[[i]] <- glmQLFTest(fit, contrast = CONTRASTS[,i]) %>%
        topTags(n=Inf) %>% as.data.frame() %>%
        rownames_to_column("Geneid")
}


bind_rows(results, .id="treatment") %>%
    write_tsv(f_out)


# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
#       mart = ensembl)

# write_tsv(Genes, "../data/Genes.list.txt")

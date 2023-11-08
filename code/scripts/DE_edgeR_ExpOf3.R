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
                 "../code/featureCounts/AllSamples_Counts.txt ../code/scratch/DE_tests.txt.gz ../code/scratch/DE_tests.mat.counts.gz", what='character')
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
protein_coding_genes <- read_tsv("/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ExpressionAnalysis/AllProteinCodingGenes.txt", col_names=c("ensembl_gene_id"))

counts <- read_tsv(f_in, comment="#") %>%
  rename_at(vars(-(1:6)),  rename_STAR_alignment_samples)

counts.ProteinCodingGenes <- counts %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  filter(Geneid %in% protein_coding_genes$ensembl_gene_id) %>%
  mutate(ensembl_gene_id = str_replace(Geneid, "(.+?)\\.\\d+$", "\\1")) %>%
  left_join(genes, by="ensembl_gene_id") %>%
  mutate(Geneid = paste(Geneid, hgnc_symbol, sep="_")) %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-ensembl_gene_id, -hgnc_symbol)

counts.ProteinCodingGenes %>%
  rownames_to_column("Geneid") %>%
  write_tsv(count_out)

Metadata <- read_tsv("../code/config/samples.3MoleculesOfInterestExperiment.tsv") %>%
  distinct(SampleName) %>%
  separate(SampleName, into=c("Experiment", "Treatment","dose.nM", "rep"), convert=T, remove=F, sep="_") %>%
  unite(Contrast, Treatment, dose.nM) %>%
  mutate(Contrast = if_else(Contrast=="DMSO_NA", "DMSO", Contrast))

counts.ProteinCodingGenes.ExpOf3 <- counts.ProteinCodingGenes %>%
  dplyr::select(Metadata$SampleName)

group <- Metadata %>% pull(Contrast) %>% as.factor() %>%
  relevel("DMSO")

design <- model.matrix(~group)

y <- DGEList(counts.ProteinCodingGenes.ExpOf3, group = group)

y$samples

CONTRASTS <- colnames(design)[-1] %>%
  makeContrasts(contrasts=., levels = design )

keep <- filterByExpr(y)
filtered <- y[keep, , keep.lib.sizes=FALSE] %>%
  calcNormFactors()
#
# filtered$samples %>%
#   rownames_to_column("sample") %>%
#   arrange(group) %>%
# ggplot(aes(x=sample, y=lib.size/1E6, fill=group)) +
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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


results.df <- bind_rows(results, .id="treatment")

write_tsv(results.df, f_out)

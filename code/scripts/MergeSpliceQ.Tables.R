#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : MergeSpliceQ.Tables
# @created     : Thursday Sep 15, 2022 13:52:20 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                # "SplicingAnalysis/SpliceQ/MergedTable. SplicingAnalysis/SpliceQ/Merged/A10-1.txt.gz SplicingAnalysis/SpliceQ/Merged/A10-2.txt.gz SplicingAnalysis/SpliceQ/Merged/A10-3.txt.gz SplicingAnalysis/SpliceQ/Merged/C2-1.txt.gz SplicingAnalysis/SpliceQ/Merged/C2-2.txt.gz SplicingAnalysis/SpliceQ/Merged/C2-3.txt.gz SplicingAnalysis/SpliceQ/Merged/CUOME-1.txt.gz SplicingAnalysis/SpliceQ/Merged/CUOME-2.txt.gz SplicingAnalysis/SpliceQ/Merged/CUOME-3.txt.gz SplicingAnalysis/SpliceQ/Merged/DMSO-1.txt.gz SplicingAnalysis/SpliceQ/Merged/DMSO-2.txt.gz SplicingAnalysis/SpliceQ/Merged/DMSO-3.txt.gz SplicingAnalysis/SpliceQ/Merged/SM2-1.txt.gz SplicingAnalysis/SpliceQ/Merged/SM2-2.txt.gz SplicingAnalysis/SpliceQ/Merged/SM2-3.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-1.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-2.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-3.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-4.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-5.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-6.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-7.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpBran-8.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-1.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-2.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-3.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-4.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-5.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-6.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-7.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpC2C5-8.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpDMSO-1.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpDMSO-2.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpDMSO-3.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-1.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-2.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-3.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-4.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-5.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-6.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-7.txt.gz SplicingAnalysis/SpliceQ/Merged/TitrationExpRis-8.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_1.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_2.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_3.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_4.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_5.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_6.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_7.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_8.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_9.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_10.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_11.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_12.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_13.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_14.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_15.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_16.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_17.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_18.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_19.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_20.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_21.txt.gz",
                 "scratch/test.out. SplicingAnalysis/SpliceQ/Merged/chRNA_21.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_8.txt.gz",
                 what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

output.fn.prefix <- args[1]
input.fn <- args[-1]

library(tidyverse)

dat.in <- lapply(setNames(input.fn, input.fn), read_tsv) %>%
    bind_rows(.id="fn")

dat.in %>%
distinct(fn, chr, IStart, IEnd, strand, transcript_ID) %>%
count(fn)

# Get most common transcript for each intron.
MostCommon <- dat.in %>%
  count(chr, IStart, IEnd, strand, transcript_ID) %>%
  group_by(chr, IStart, IEnd, strand) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  dplyr::select(-n)

dat.wide <- dat.in %>%
  inner_join(MostCommon) %>%
  unite(Intron, chr, IStart, IEnd, strand) %>%
  mutate(SJCovSplit = sj5_cov_split + sj3_cov_split,
         SJCovNonSplit = sj5_cov_nonsplit+sj5_cov_split+sj3_cov_nonsplit+sj3_cov_split) %>%
  mutate(SE = SJCovSplit/SJCovNonSplit) %>%
  dplyr::select(-SJCovSplit, -SJCovNonSplit) %>%
  group_by(fn, Intron) %>%
  summarise_at(vars(exon5_cov:SE), median, na.rm = TRUE) %>%
  ungroup() %>%
  pivot_wider(names_from = "fn", values_from=exon5_cov:SE, id_cols="Intron", id_expand=T)

for (col in c(colnames(dat.in)[9:16], "SE")){
  output.fn <- paste0(output.fn.prefix, col, ".txt.gz")
  print(paste("writing", output.fn))
  dat.wide %>%
    dplyr::select(Intron, starts_with(paste0(col, "_"))) %>%
    write_tsv(output.fn)
}


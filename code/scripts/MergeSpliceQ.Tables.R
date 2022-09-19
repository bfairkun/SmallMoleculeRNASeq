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
                 "scratch/test.out.txt SplicingAnalysis/SpliceQ/Merged/chRNA_9.txt.gz SplicingAnalysis/SpliceQ/Merged/chRNA_8.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

output.fn <- args[1]
input.fn <- args[-1]

library(tidyverse)

dat.in <- lapply(setNames(input.fn, input.fn), read_tsv) %>%
    bind_rows(.id="fn")

dat.in %>%
distinct(fn, chr, IStart, IEnd, strand, transcript_ID) %>%
count(fn)

dat.wide <- dat.in %>%
  unite(Intron, chr, IStart, IEnd, strand) %>%
  mutate(SE.perrow = (sj5_cov_split + sj3_cov_split)/(sj5_cov_nonsplit+sj5_cov_split+sj3_cov_nonsplit+sj3_cov_split)) %>%
  group_by(fn, Intron) %>%
  summarise(SE = median(SE.perrow),
            IER=median(IER)) %>%
  ungroup() %>%
  pivot_wider(names_from = "fn", values_from=c("SE", "IER"), id_cols="Intron", id_expand=T)

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
    distinct(fn, chr, IStart, IEnd, strand, transcript_ID)
    # count(fn)

dat.in %>%
    distinct(fn, chr, IStart, IEnd, strand, exon5_cov, exon3_cov)

dat.in %>% 
    head() %>%
    rowwise() %>% 
    mutate(sumrange = sum(c_across(exon5_cov:exon3_cov), na.rm = T)) %>%
    filter(sumrange >= 10)

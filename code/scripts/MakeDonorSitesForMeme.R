#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : MakeDonorSitesForMeme
# @created     : Wednesday Sep 27, 2023 10:01:52 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../output/EC50Estimtes.FromPSI.txt.gz DonorMotifSearches/Sites/", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


f_in <- args[1]
f_out_prefix <- args[2]

library(tidyverse)

dat <- read_tsv(f_in)

BranaplamAndRisdiplamDonorSequences <- dat %>%
  filter(!is.na(Steepness)) %>%
  filter(EC.Ratio.Test.Estimate.FDR_Branaplam.Risdiplam < 0.1) %>%
  mutate(BranaplamOrRisdiplamSite = if_else(ECRatio.ComparedToGenomewideMedian_Branaplam.Risdiplam > 1, "RisdiplamSpecific", "BranaplamSpecific")) %>%
  # count(BranaplamOrRisdiplamSite) %>%
  dplyr::select(BranaplamOrRisdiplamSite, seq, junc)
  
BranaplamAndRisdiplamDonorSequences %>%
    dplyr::select(seq, BranaplamOrRisdiplamSite) %>%
    group_by(BranaplamOrRisdiplamSite) %>%
   group_walk(~ write_tsv(.x, paste0(f_out_prefix, .y$BranaplamOrRisdiplamSite, ".txt"), col_names=F))

#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : GetBranchWindow
# @created     : Friday Feb 11, 2022 15:29:09 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../data/BradleyDistTo3ss.txt.gz scratch/memetest.50.fa.tab Meme/PSP/BP.psp Meme/Fasta/BP.fa 5", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

bradley_bp_dist_in <- args[1]
f_in <- args[2]
f_psp_out <- args[3]
f_fa_out <- args[4]
w <- as.numeric(args[5])

library(tidyverse)

bp_dist <- read_delim(bradley_bp_dist_in, delim=' ', col_names = c("intron_type", "dist"))
dat <- read_tsv(f_in, col_names=c("Name", "Seq"))
write_lines(c(">BP 5"), f_out)

priors <- 
bp_dist %>%
    filter(intron_type == "u2") %>%
    # filter(dist < 100 & dist > 5) %>%
    filter(dist <= 100 & dist >= 8) %>%
    count(dist) %>%
    mutate(n= lag(n, 2)) %>%
    right_join(data.frame(dist=1:50), by="dist") %>%
    replace_na(list(n=0)) %>%
    arrange(desc(dist)) %>%
    mutate(rn = row_number()) %>%
    mutate(n = rev(n)) %>%
    mutate(n = case_when(
                         rn < w ~ 0,
                         TRUE ~ n
                         )) %>%
    mutate(prior = n/sum(n)) %>%
    # ggplot(aes(x=dist, y=prior)) +
    # geom_line()
    arrange(dist) %>% pull(prior)
priors

dat %>%
    mutate(str_out = str_glue(">{Name}\n{paste(priors[1:45], collapse=' ')}")) %>%
    select(str_out) %>%
    write.table(f_psp_out, col.names=F, quote=F, row.names=F)

dat %>%
    mutate(Seq = substr(Seq, 0, 45)) %>%
    mutate(str_out = str_glue(">{Name}\n{Seq}")) %>%
    select(str_out) %>%
    write.table(f_fa_out, col.names=F, quote=F, row.names=F)

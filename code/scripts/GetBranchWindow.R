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
                 "../data/BradleyDistTo3ss.txt.gz Meme/PSP/BP.txt", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_in <- args[1]
f_out <- args[2]

library(tidyverse)

bp_dist <- read_delim(f_in, delim=' ', col_names = c("intron_type", "dist"))

bp_dist %>%
    filter(intron_type == "u2") %>%
    filter(dist < 100 & dist > 5) %>%
    pull(dist) %>%
    hist()
    # quantile(c(0.05, 0.95))

bp_dist %>%
    filter(intron_type == "u2") %>%
    # filter(dist < 100 & dist > 5) %>%
    filter(dist < 95 & dist > 5) %>%
    ggplot(aes(x=dist)) +
    geom_density(bw=3)
    # quantile(c(0.05, 0.95))

bp_dist %>%
    filter(intron_type == "u2") %>%
    # filter(dist < 100 & dist > 5) %>%
    filter(dist <= 50 & dist >= 5) %>%
    ggplot(aes(x=dist)) +
    geom_bar()
    # quantile(c(0.05, 0.95))


write_lines(c(">BP 5"), f_out)

bp_dist %>%
    filter(intron_type == "u2") %>%
    # filter(dist < 100 & dist > 5) %>%
    filter(dist <= 50 & dist >= 5) %>%
    count(dist) %>%
    full_join(data.frame(dist=1:50), by="dist") %>%
    arrange(dist) %>%
    replace_na(list(n=0)) %>%
    pull(n) %>% paste(collapse = ' ') %>%
    write_lines(f_out, append=T)
    

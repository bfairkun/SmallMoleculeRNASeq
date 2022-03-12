
# input bed file should already be annotated by repeats and exon/intron
# this script just classify repeats to repeat families


library(tidyverse)

args = commandArgs(trailingOnly = T)
lookup = read_tsv(args[2])
#bed = data.table::fread("/project2/yangili1/cdai/rna_edit/chrRNA/Results/hs38/BCFCall/FinalOutput/test.bed")
bed = data.table::fread(args[1])

left_join(bed, lookup[,c(1,3)], by = c("V11" = "repName")) %>%
  mutate_at("repFamily", ~ replace_na(.x, ".")) %>%
  write_tsv(args[1], col_names = F)

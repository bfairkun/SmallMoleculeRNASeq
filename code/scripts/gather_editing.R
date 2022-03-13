library(tidyverse)
library(data.table)

runmode = 1

if (runmode == 1) {
    input = snakemake@input
    output.bed = snakemake@output[[1]]
    ap_min = snakemake@params[['ap_min']]
    N_min = snakemake@params[['N_min']]
} else {
    setwd("/project2/yangili1/cdai/mRNA-editing/YRI")
    input = dir('Results/hs38/EditLevel', '^NA[0-9]+\\.txt$', full.names = T)
    ap_min = 1
    N_min = 10
}

names(input) = str_split(input, '/', simplify = T)[, 3] %>% str_remove("\\.txt")

# read in files
edits = map(input, ~ fread(.x, sep = "\t", header = T))
# remove NA in EL col
edits = map(edits, ~ mutate_at(.x, "EL", function(v) replace_na(v, 0)))

# make count matrix of editing site read counts
AP_mx = imap(edits, ~add_column(.x, SAMP = .y) %>% unique) %>% 
    bind_rows() %>% 
    select(chr:BEDend, SAMP, AP) %>% 
    pivot_wider(names_from = c('SAMP'), values_from = c('AP'), values_fill = 0) %>% 
    unique %>% 
    mutate(site = paste(chr, BEDstart, sep = ":"), .after = "BEDend") %>% 
    select(-chr, -BEDstart, -BEDend) %>% 
    column_to_rownames("site")

# make editing level matrix
EL_mx = imap(edits, ~add_column(.x, SAMP = .y) %>% unique) %>% 
    bind_rows() %>% 
    select(chr:BEDend, SAMP, EL) %>% 
    pivot_wider(names_from = c('SAMP'), values_from = c('EL'), values_fill = 0) %>% 
    unique %>% 
    mutate(site = paste(chr, BEDstart, sep = ":"), .after = "BEDend") %>% 
    select(-chr, -BEDstart, -BEDend) %>% 
    column_to_rownames("site")


# get number of samples per min AP
getNSamplePerMinAP = function(mx, MIN) {
    mx = modify(mx, ~ .x >= MIN)  %>% 
            rowSums
    return(mx)
}

# get row index (editing sites) that passed minimum thresholds
keep_idx = getNSamplePerMinAP(AP_mx, ap_min)
keep_idx = which(keep_idx >= N_min)
sample_names = colnames(EL_mx)

# selected sites that passed
sites = rownames(EL_mx[keep_idx, ]) %>% 
    str_split(":", simplify = T)

# select filtered sites and make phenotype bed
out_df = EL_mx[keep_idx, ] %>%
    rownames_to_column("pid") %>% # make QTL required Bed format
    mutate(chr = sites[,1],
           BEDstart = as.integer(sites[,2]),
           BEDend = as.integer(sites[,2]) + 1,
           .before = "pid") %>% 
    mutate(gid = pid, strand = ".", .after = "pid")

# sort chromosome and start
chroms = paste("chr", c(as.character(1:22), "X", "Y"), sep = "")
chroms = factor(chroms, chroms)
out_df = mutate_at(out_df, "chr", ~factor(.x, levels = chroms)) %>%
    arrange(chr, BEDstart, BEDend) %>%
    mutate_at(c("BEDstart", "BEDend"), as.integer)

# fix bed header
colnames(out_df) = c("#Chr", colnames(out_df)[-1])

# write output 
fwrite(out_df, output.bed, sep = "\t", col.names = T, row.names = F)

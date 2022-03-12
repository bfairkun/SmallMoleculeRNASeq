library(tidyverse)

runmode = 1

if (runmode == 1) {
    input_counts = snakemake@input[['cnt']]
    input_beds = snakemake@input[['bed']]
    output = snakemake@output[[1]]
} else {
    setwd("/project2/yangili1/cdai/mRNA-editing/YRI")
    input_counts = "Results/hs38/countCoverage/NA18486.txt"
    input_beds = "Results/hs38/BCFCall/raw/NA18486.bed"
    output = "test.txt"
}


counts = data.table::fread(input_counts, sep = "\t", header = TRUE)

# read bed that includes mismatch type
beds = data.table::fread(input_beds, sep = "\t", header = FALSE)
# keep only A>G (T>C) mismatches
beds = filter(beds, V4 %in% c("T,C", "A,G"))

# complement counts with mismatch type
counts = left_join(counts, beds[, 1:4],
            by = c("chr" = "V1", "BEDstart" = "V2", "BEDend" = "V3"))  %>%
            drop_na

# copute DP, AP, and EL (Editing Level)

counts = rename(counts, EditType = V4)  %>%
            mutate(DP = A + C + G + T,
                   AP = case_when(EditType == "A,G" ~ G,
                                  EditType == "T,C" ~ C))  %>% 
            mutate(EL = round(AP/DP, 5))  %>%
            select(-A, -C, -G, -T) %>%
            mutate_at(c("BEDstart", "BEDend"), as.integer)

# write to tsv
write_tsv(counts, output)
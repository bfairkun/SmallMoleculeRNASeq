library(tidyverse)


# set run mode: 1 - snakemake ; 2 - debug in rstudio directly
runmode = 1

# inputfile is tab delimited output in 
if (runmode == 1) {
  input.file = snakemake@input[[1]]
  output.file = snakemake@output[[1]]
} else {
  input.file = "Results/hs38/VariantsBCF/checkDB/ERR188040.bed"
  output.file = "ERR188040_test.bed"
}


  COL_names = c("chr", "start", "end", "id", "score", "strand", 
              "ref", "alt", "qual", "dp", "ad_ref", "ad_alt",
              "ef", "rnadb", "mapq")

BED = data.table::fread(input.file, header = F, col.names = COL_names)


# hard criteria
MIN_MISMATCH = 1
#AF_THREASHOLD = 0.5

BED_filtered = filter(BED, ad_alt >= MIN_MISMATCH)

BED_filtered = mutate(BED_filtered, 
                      label = if_else(str_detect(rnadb, "RADAR"),
                                      "known", "novel"))

select(BED_filtered, chr:end, label, ef, strand, ad_ref, ad_alt) %>% 
  mutate_at("ef", ~ round(.x*1000)) %>% 
  data.table::fwrite(output.file, sep = "\t", quote = F, col.names = F)

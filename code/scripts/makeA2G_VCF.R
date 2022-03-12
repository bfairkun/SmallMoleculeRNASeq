##### Gather sites
library(tidyverse)


if (exists("snakemake")) {
  print("### taking snakemake inputs")
  input.files.ls = snakemake@input # note this is a list
  output.file.name = snakemake@output[[1]]
} else {
  print("### working in Debug mode, not taking snakemake inputs!!!")
  input.files.ls = dir("Results/hs38/BCFCall/FinalAnno", "^JS.+[0-9]+\\.bed", full.names = T)
  output.file.name = "Results/hs38/BCFCall/FinalOutput/STAR-WASP-A2G-editing-remap.vcf"
}


# Load data files ---------------------------------------------------------

chBED.ls = input.files.ls
names(chBED.ls) = chBED.ls %>%
  str_split("/", simplify = T) %>%
  .[, ncol(.)] %>%
  str_replace_all("(JS\\-YZ\\-10S\\-)|(\\.bed)|(\\-)", "")
COLS = c("chr", "start", "end", "label", "af", "strand",
        "ref_ad", "alt_ad", "gene", "feature", "repFam", "repCla")
chBED.ls = map(chBED.ls, ~ data.table::fread(
            .x, sep = "\t", header = F, col.names = COLS))




# Union all A>G sites ---------------------------------------------------

# filter editing sites
# effectively not filtering anything
af.cutoff = 0
alt_ad.cutoff = 0

chBED.filtered.ls = map(chBED.ls,
    ~ filter(.x, af > af.cutoff & alt_ad > alt_ad.cutoff))

unique_sites = imap(chBED.filtered.ls,
    function(df, df.nm) select(df, chr, start, end)) %>%
  do.call(bind_rows, .) %>%
  unique

# order by chr then by start coordinates
# a bit clunky because arrange does not do natural sort
unique_sites = group_by(unique_sites, chr) %>%
                  group_split %>% 
                  map_df(~arrange(.x, start)) %>%
                  .[gtools::mixedorder(.$chr), ]


# write output ------------------------------------------------------------


vcf.out = data.frame(CHROM=unique_sites$chr, POS=unique_sites$start,
                     ID=".", REF="A", ALT="G", QUAL=".", FILTER="PASS",
                     INFO=".", GT="GT", "SPL"="0/1")

colnames(vcf.out) = c("#CHROM","POS",'ID','REF','ALT','QUAL','FILTER','INFO','GT','SPL')

data.table::fwrite(vcf.out, file = output.file.name, quote = F, sep = "\t")

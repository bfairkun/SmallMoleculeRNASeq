require(tidyverse)
require(GenomicRanges)

if (interactive()) {
    args = scan(text="Data/BF_Alignment/NA19108/1/SJ.out.tab  Data/BF_Alignment/NA19108/2/SJ.out.tab",
                what = "character")
} else {
    args = commandArgs(trailingOnly = TRUE)
}


FileIn = args[1:length(args)-1]
FileOut = args[length(args)]

print(paste("Input: ", FileIn))
print(paste("Output: ", FileOut))

# refer to STAR manual for sj output file columns
sj.cols = c('chr','intron_st','intron_en','strand','sj_motif',
            'annot_flag','unique_reads','mul_reads','max_sj_over')
sj.tab = map_dfr(FileIn, ~read_tsv(.x, col_names = sj.cols, col_types = "ciiiiiiii"))

# combine read counts from all ERR files for each sample
sj.tab = group_by(sj.tab, chr, intron_st, intron_en, strand, sj_motif, annot_flag) %>% 
    summarise(unique_reads = sum(unique_reads),
              mul_reads = sum(mul_reads),
              max_sj_over = sum(max_sj_over), .groups = "drop")

# set 20 as hard filter for combined unique reads and multi-mapping reads span SJ
MIN_READS = 20
sj.tab = mutate(sj.tab, com_reads = unique_reads + mul_reads) %>% 
            filter(com_reads >= MIN_READS & str_detect(chr, "chr[0-9XY][0-9]?"))

sj.tab = mutate(sj.tab, strand2 = case_when(
    strand == 0 ~ ".",
    strand == 1 ~ "+", 
    strand == 2 ~ "-"
))



# get the coordinates of 4bp intron region that is: 
# 5'***---------intron-----------***3' produces
# coordinates for 5'***
# and coordinates for 3'***

# convert STAR produced SJ sites into genomic ranges
sj.grange = makeGRangesFromDataFrame(sj.tab, keep.extra.columns = F, ignore.strand = F, 
                                    seqnames.field = 'chr', start.field = 'intron_st', 
                                    end.field = 'intron_en',
                                    strand.field = "strand2")

# 4bp in introns immediately adjacent to the 5' splice site
sj.5p = resize(sj.grange, width = 4, fix = "start")

# 4bp in introns immediately adjacent to the 3' splice site
sj.3p = resize(sj.grange, width = 4, fix = "end")

# covnert to BED format 0 based left closed, right open
sj = sort(c(sj.5p, sj.3p)) %>% as.data.frame %>% 
    mutate_at('start', ~ .x - 1) %>% 
    mutate_at('strand', ~ str_replace(.x, "\\*", "\\.")) %>% 
    add_column(name=".", score=".") %>% 
    select(seqnames, start, end, name, score, strand)

print(paste("# Saving combined SJ 4bp BED files to: ", FileOut))
write_tsv(sj, FileOut, col_names = F)

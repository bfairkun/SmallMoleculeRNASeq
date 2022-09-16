#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TidyAndSpearmanCorTitrationData
# @created     : Thursday Sep 15, 2022 18:48:24 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(edgeR)

sample.list <- read_tsv("bigwigs/BigwigList.tsv",
                        col_names = c("SampleName", "bigwig", "group", "strand")) %>%
  filter(strand==".") %>%
  mutate(old.sample.name = str_replace(bigwig, "/project2/yangili1/bjf79/20211209_JingxinRNAseq/code/bigwigs/unstranded/(.+?).bw", "\\1")) %>%
  separate(SampleName, into=c("treatment", "dose.nM", "cell.type", "libType", "rep"), convert=T, remove=F, sep="_") %>%
  left_join(
    read_tsv("bigwigs/BigwigList.groups.tsv", col_names = c("group", "color", "bed", "supergroup")),
    by="group"
  ) %>%
  dplyr::select(-strand, -bed, -bigwig)

all.samples.PSI <- read_tsv("SplicingAnalysis/leafcutter_all_samples/PSI.table.bed.gz")

titration.series.counts <- read_tsv("featureCounts/Counts.titration_series.txt", comment="#") %>%
  rename_at(vars(-(1:6)), ~str_replace(.x, "Alignments/STAR_Align/(.+?)/Aligned.sortedByCoord.out.bam", "\\1"))

ExpressedGenes.cpm <- titration.series.counts %>%
  dplyr::select(Geneid, everything(), -c(2:6)) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  cpm(prior.count=0.1) %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    gather(key="old.sample.name", value="CPM", -Geneid) %>%
    inner_join(sample.list, by="old.sample.name")

all.samples.5ss <- read_tsv("SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz", col_names = c("intron", "seq", "score")) %>%
  mutate(intron = str_replace(intron, "^(.+?)::.+$", "\\1")) %>%
  separate(intron, into=c("chrom", "start", "stop", "strand"), sep="_", convert=T, remove=F)

all.samples.intron.annotations <- read_tsv("SplicingAnalysis/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz") %>%
  mutate(IntronType = recode(anchor, A="Alt 5'ss", D="Alt 3'ss", DA="Annotated", N="New Intron", NDA="New Splice Site combination"))


merged.dat.df <- all.samples.PSI %>%
  dplyr::select(1:6, contains("LCL")) %>%
  dplyr::select(-contains("chRNA")) %>%
  drop_na() %>%
  inner_join(all.samples.5ss, by=c("#Chrom"="chrom", "start", "end"="stop", "strand")) %>%
  inner_join(all.samples.intron.annotations, by=c("#Chrom"="chrom", "start", "end", "strand") )

set.seed(0)
dmso.data <- merged.dat.df %>%
  dplyr::select(junc, contains("LCL")) %>%
  gather("SampleName", "PSI", -junc) %>%
  left_join(sample.list, by="SampleName") %>%
  filter(treatment=="DMSO") %>%
  arrange(junc, rep) %>%
  mutate(dose.nM = 0)
dmso.data$treatment <- replicate(nrow(dmso.data)/3, sample(c("Branaplam", "Risdiplam", "C2C5"), 3)) %>% as.vector()

# Tidy the splicing.gagt.df dataframe
merged.dat.df.tidy <- merged.dat.df %>%
  dplyr::select(junc, contains("LCL")) %>%
  gather("SampleName", "PSI", contains("LCL")) %>%
  left_join(sample.list, by="SampleName") %>%
  filter(!treatment=="DMSO") %>%
  bind_rows(dmso.data) %>%
  inner_join(
    merged.dat.df %>%
        dplyr::select(1:6, seq, Donor.score=score.x, gene_names, gene_ids, IntronType),
    by="junc"
  )


spearman.coefs <- merged.dat.df.tidy %>%
    group_by(treatment, junc) %>%
    summarise(spearman = cor(dose.nM, PSI, method='sp'))

# Write out splicing results
merged.dat.df.tidy %>%
    inner_join(spearman.coefs, by=c("treatment", "junc")) %>%
    write_tsv("DoseResponseData/LCL/TidySplicingDoseData.txt.gz")
  

# Tidy the cpm df
set.seed(0)
dmso.data <- ExpressedGenes.cpm %>%
  filter(treatment=="DMSO") %>%
  arrange(Geneid, rep) %>%
  mutate(dose.nM = 0)
dmso.data$treatment <- replicate(nrow(dmso.data)/3, sample(c("Branaplam", "Risdiplam", "C2C5"), 3)) %>% as.vector()

ExpressedGenes.cpm.tidy <-
    ExpressedGenes.cpm %>%
    filter(!treatment=="DMSO") %>%
    bind_rows(dmso.data)

spearman.coefs <- ExpressedGenes.cpm.tidy %>%
    group_by(treatment, Geneid) %>%
    summarise(spearman = cor(dose.nM, CPM, method='sp'))

#Write out expression results
ExpressedGenes.cpm.tidy %>%
    inner_join(spearman.coefs, by=c("treatment", "Geneid")) %>%
    write_tsv("DoseResponseData/LCL/TidyExpressionDoseData.txt.gz")

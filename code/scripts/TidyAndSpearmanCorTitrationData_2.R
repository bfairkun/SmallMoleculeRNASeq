#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TidyAndSpearmanCorTitrationData
# @created     : Thursday Sep 15, 2022 18:48:24 CDT
#
# @description : Use previous tidied data and keep all 3 DMSO replicates for each treatment series. Divide PSI by within-cluster within-sample max PSI for easier interpretability... And also convert CPM to log2FC relative to DMSO.
# These are the metrics I will use for modelling.
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
library(data.table)

splicing.dat <- fread("DoseResponseData/LCL/TidySplicingDoseData.txt.gz") %>%
  as_tibble()

splicing.dat.dmso <- splicing.dat %>%
  filter(dose.nM == 0)

splicing.dat.dmso.Replicated <- bind_rows(
  splicing.dat.dmso %>% mutate(treatment="C2C5"),
  splicing.dat.dmso %>% mutate(treatment="Risdiplam"),
  splicing.dat.dmso %>% mutate(treatment="Branaplam")
)

splicing.dat.dmso.Replicated %>%
  dplyr::select(junc, gid)

splicing.dat.Modified <- bind_rows(
  splicing.dat %>% filter(!dose.nM == 0),
  splicing.dat.dmso.Replicated
) %>%
  group_by(SampleName, gid) %>%
  mutate(MaxWithinSampleWithinCluster = max(PSI,na.rm=T)) %>%
  mutate(PSI_transformed = PSI/MaxWithinSampleWithinCluster) %>%
  ungroup() %>%
  group_by() %>%
  dplyr::select(-spearman)


spearman.coefs <- splicing.dat.Modified %>%
  group_by(treatment, junc) %>%
  summarise(spearman = cor(dose.nM, PSI_transformed, method='sp'))

# Write out splicing results
splicing.dat.Modified %>%
  inner_join(spearman.coefs, by=c("treatment", "junc")) %>%
  dplyr::select(-MaxWithinSampleWithinCluster) %>%
  # distinct(treatment, junc, .keep_all=T) %>%
  # group_by(treatment) %>%
  # sample_n(5000) %>%
  # ungroup() %>%
  # # dplyr::select(spearman.x, spearman.y) %>%
  # ggplot(aes(x=spearman.x, y=spearman.y)) +
  # geom_point(alpha=0.01) +
  # facet_wrap(~treatment)
  write_tsv("DoseResponseData/LCL/TidySplicingDoseData_PSITransformedAndAllDMSORepsInEachSeries.txt.gz")


# When I do modelling, I think it will make sense to simultaneously fit all treatments, and constrain limits to be the same in all three treatments... It's just ED50 and slope that can vary b/n treatments
Expression.dat <- fread("DoseResponseData/LCL/TidyExpressionDoseData.txt.gz") %>%
  as_tibble()

Expression.dat.dmso <- Expression.dat %>%
  filter(dose.nM == 0)

Expression.dat.dmso.Replicated <- bind_rows(
  Expression.dat.dmso %>% mutate(treatment="C2C5"),
  Expression.dat.dmso %>% mutate(treatment="Risdiplam"),
  Expression.dat.dmso %>% mutate(treatment="Branaplam")
)

MinCPM <- Expression.dat.dmso.Replicated %>%
  filter(CPM>0) %>%
  pull(CPM) %>% min()

MinCPM

MeanDMSO.Expression <- Expression.dat.dmso.Replicated %>%
  group_by(Geneid) %>%
  filter(all(CPM>0)) %>%
  ungroup() %>%
  # filter out genes with any DMSO samples with 0 counts
  #obviates need for psueodcounts, but perhaps too stringent?
  #... 21,167 genes remaining, so probably not too stringent.
  # distinct(Geneid)
  mutate(CPM = CPM + MinCPM/10) %>%
  group_by(Geneid) %>%
  summarise(Log2CPM_Mean_DMSO = mean(log2(CPM), na.rm=T)) %>%
  ungroup()

Symbols <- read_tsv("../data/Genes.list.txt")

Expression.dat.Modified <- bind_rows(
  Expression.dat %>% filter(!dose.nM == 0),
  Expression.dat.dmso.Replicated
) %>%
  left_join(MeanDMSO.Expression) %>%
  mutate(log2CPM = log2(CPM+MinCPM/10)) %>%
  mutate(log2FC=log2CPM-Log2CPM_Mean_DMSO) %>%
  mutate(ensembl_gene_id = str_replace(Geneid, "(^E.+?)\\..+$", "\\1")) %>%
  # dplyr::select(Geneid, ensembl_gene_id)
  left_join(Symbols) %>%
  dplyr::select(-spearman)



Expression.dat.Modified %>%
  filter(hgnc_symbol %in% c("SMN2", "HTT", "STAT1", "GALC", "ATP5", "SMN1")) %>%
  ggplot(aes(x=dose.nM, y=log2FC, color=treatment)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(trans='log1p') +
  facet_wrap(~hgnc_symbol)

spearman.coefs <- Expression.dat.Modified %>%
  group_by(treatment, Geneid) %>%
  summarise(spearman = cor(dose.nM, log2FC, method='sp'))

Expression.dat.Modified %>%
  inner_join(spearman.coefs, by=c("treatment", "Geneid")) %>%
  write_tsv("DoseResponseData/LCL/TidyExpressionDoseData_logFCTransformedAndAllDMSORepsInEachSeries.txt.gz")





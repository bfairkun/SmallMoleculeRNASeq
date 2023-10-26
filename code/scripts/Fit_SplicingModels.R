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
                 "DoseResponseData/LCL/TidySplicingDoseData_PSITransformedAndAllDMSORepsInEachSeries.txt.gz scratch/FitSplicingModels.RData", what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}

f_in <- args[1]
f_out <- args[2]

library(tidyverse)
library(data.table)
library(drc)

sample_n_of <- function(data, size, ...) {
  dots <- quos(...)

  group_ids <- data %>%
    group_by(!!! dots) %>%
    group_indices()

  sampled_groups <- sample(unique(group_ids), size)

  data %>%
    filter(group_ids %in% sampled_groups)
}

splicing.dat <- fread(f_in) %>%
  group_by(treatment) %>%
  mutate(doseRank = dense_rank(dose.nM)) %>%
  ungroup() %>%
  as_tibble() %>%
  mutate(doseInRisdiscale = case_when(
    treatment == "C2C5" ~ dose.nM * 10,
    treatment == "Branaplam" ~ dose.nM * sqrt(10),
    TRUE ~ dose.nM
  ))

splicing.dat %>%
  distinct(doseRank, dose.nM, treatment) %>%
  mutate(Log10dose.nM=log10(dose.nM)) %>%
  pivot_wider(names_from="treatment", values_from=c("Log10dose.nM", "dose.nM"))

splicing.dat %>%
  filter(!is.na(spearman) & !is.na(PSI_transformed)) %>%
  add_count(junc) %>%
  filter(n==33) %>%
  group_by(junc) %>%
  filter(any(abs(spearman)>0.5)) %>%
  ungroup() %>%
  sample_n_of(40, junc) %>%
  ggplot(aes(x=doseRank, y=PSI_transformed, color=treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~junc, scale="free_y")



model.dat.df.PSI <-
  splicing.dat %>%
  filter(!is.na(spearman) & !is.na(PSI_transformed)) %>%
  add_count(junc) %>%
  filter(n==33) %>%
  group_by(junc) %>%
  filter(any(abs(spearman)>0.5)) %>%
  ungroup() %>%
  nest(data=-c("junc")) %>%
  sample_n_of(40, junc)

model.dat.df.PSI$junc[1] %>% as.data.frame()

Results.models <- list()
for(i in 1:nrow(model.dat.df.PSI)) {
  # for(i in 1:100) {
  tryCatch(
    expr = {
      Geneid <- model.dat.df.PSI$junc[i]
      data <- model.dat.df.PSI$data[i] %>% as.data.frame()
      fit <- drm(formula = PSI_transformed ~ doseInRisdiscale,
                 data = data,
                 fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,NA,NA,NA)),
                 lowerl = c(-Inf, 0, -Inf, -Inf),
                 upperl = c(Inf, 100, Inf, Inf),
                 # fct = LL.4(),
                 curveid = treatment,
                 # pmodels=data.frame(1, 1, treatment, treatment),
                 pmodels=data.frame(treatment, 1, treatment, treatment),
                 robust = "mean"
      )
      Results.models[[Geneid]] <- fit
      message("Successfully fitted model.")
    },
    error=function(e){
      if (i < 100){
        cat("ERROR :",conditionMessage(e), Geneid, "\n")
      }
    })
}


#Examples that failed to converge
splicing.dat %>%
  filter(junc %in% c("chr8:73009205:73020684:clu_17520_+", "chr5:138444867:138445234:clu_11395_+", "chr5:88722829:88728493:clu_12024_-", "chr5:79306461:79312403:clu_11118_+", "chr4:1179541:1201269:clu_10055_-", "chr3:58433684:58433768:clu_8528_-", "chr22:40964117:40967799:clu_38775_+", "chr2:210017559:210018153:clu_6621_+", "chr2:98610895:98616432:clu_6039_+", "chr19:4055247:4067141:clu_34156_-")) %>%
  ggplot(aes(x=doseInRisdiscale, y=PSI_transformed, color=treatment)) +
  geom_line() +
  scale_x_continuous(trans='log1p') +
  facet_wrap(~junc, scale="free_y")

#Examples that did converge
splicing.dat %>%
  filter(junc %in% (Results.models %>% names())) %>%
  ggplot(aes(x=doseInRisdiscale, y=PSI_transformed, color=treatment)) +
  geom_line() +
  scale_x_continuous(trans='log1p') +
  facet_wrap(~junc, scale="free_y")

Results.models$`chr9:77382087:77384585:clu_18943_+` %>% plot()

Results.models$`chr20:528049:543672:clu_36498_-` %>% plot()

save(Results.models, file = f_out)

# plot(Results.models$ENSG00000135148.12)
# summary(Results.models$ENSG00000135148.12)
#
# plot(Results.models$ENSG00000115415.20)
# summary(Results.models$ENSG00000115415.20)


# getData <- function(fit){
#   return(fit$data)
# }
#
# dat <- lapply(Results.models, getData ) %>%
#   bind_rows(.id="gene")


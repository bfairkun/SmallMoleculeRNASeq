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
                 "DoseResponseData/LCL/TidyExpressionDoseData_logFCTransformedAndAllDMSORepsInEachSeries.txt.gz scratch/FitExpressionModels.RData", what='character')
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

expression.dat <- fread(f_in) %>%
  group_by(treatment) %>%
  mutate(doseRank = dense_rank(dose.nM)) %>%
  ungroup() %>%
  as_tibble() %>%
  mutate(doseInRisdiscale = case_when(
    treatment == "C2C5" ~ dose.nM * 10,
    treatment == "Branaplam" ~ dose.nM * sqrt(10),
    TRUE ~ dose.nM
  ))

expression.dat %>%
  distinct(doseRank, dose.nM, treatment) %>%
  mutate(Log10dose.nM=log10(dose.nM)) %>%
  pivot_wider(names_from="treatment", values_from=c("Log10dose.nM", "dose.nM"))

expression.dat %>%
  filter(!is.na(spearman) & !is.na(log2FC)) %>%
  add_count(Geneid) %>%
  filter(n==33) %>%
  sample_n_of(40, hgnc_symbol) %>%
  ggplot(aes(x=doseRank, y=log2FC, color=treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~hgnc_symbol, scale="free_y")

GenesToHighlight <- c("STAT1", "HTT", "MYB",  "TRIM11", "TRAFD1", "VEGFA", "FBXW11", "HSD17B4")

expression.dat %>%
  mutate(treatment = factor(treatment)) %>%
  filter(hgnc_symbol %in% GenesToHighlight) %>%
  ggplot(aes(x=doseRank, y=log2FC, color=treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~hgnc_symbol, scale="free_y")


expression.dat %>%
  mutate(treatment = factor(treatment)) %>%
  filter(hgnc_symbol %in% GenesToHighlight) %>%
  ggplot(aes(x=doseRank, y=log2FC, color=treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~hgnc_symbol)

model.dat.df.CPM.StandardNormFactors <-
  expression.dat %>%
  filter(!is.na(spearman) & !is.na(log2FC)) %>%
  add_count(Geneid) %>%
  filter(n==33) %>%
  mutate(treatment = factor(treatment)) %>%
  filter(hgnc_symbol %in% GenesToHighlight) %>%
  nest(data=-c("Geneid"))

Results.models <- list()
for(i in 1:nrow(model.dat.df.CPM.StandardNormFactors)) {
  # for(i in 1:100) {
  tryCatch(
    expr = {
      Geneid <- model.dat.df.CPM.StandardNormFactors$Geneid[i]
      data <- model.dat.df.CPM.StandardNormFactors$data[i] %>% as.data.frame()
      SignSpearman <- (sign(median(sign(data$spearman))))
      if (SignSpearman == 1){
        FixedParams <- c(NA,0,NA,NA)
        ParamNames <- c("Steepness", "LimitAtDoseZero", "LimitAtDoseInfinity", "ED50")
        # pmodels <- data.frame(treatment, treatment, 1, treatment)
      }
      else{
        FixedParams <- c(NA,NA,0,NA)
        ParamNames <- c("Steepness", "LimitAtDoseInfinity", "LimitAtDoseZero", "ED50")
        # pmodels <- data.frame(treatment, treatment, 1, treatment)
      }
      fit <- drm(formula = log2FC ~ doseInRisdiscale,
                 data = data,
                 fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,NA,0,NA)),
                 # fct = LL.4(names=ParamNames, fixed=FixedParams),
                 # fct = LL.4(),
                 curveid = treatment,
                 # pmodels=data.frame(1, 1, treatment, treatment),
                 # pmodels=data.frame(treatment, 1, treatment, treatment),
                 pmodels=data.frame(treatment, 1, 1, treatment),
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

# save(Results.models, file = f_out)


plot(Results.models$ENSG00000135148.12)
summary(Results.models$ENSG00000135148.12)

plot(Results.models$ENSG00000115415.20)
summary(Results.models$ENSG00000115415.20)


expression.dat %>%
  filter(hgnc_symbol == "TRIM11") %>%
  pull(Geneid) %>% unique()

# getData <- function(fit){
#   return(fit$data)
# }
#
# dat <- lapply(Results.models, getData ) %>%
#   bind_rows(.id="gene")

expression.dat %>%
  filter(Geneid %in% c("ENSG00000154370.16", "ENSG00000115415.20")) %>%
  ggplot(aes(x=doseRank, y=log2FC, color=treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~hgnc_symbol, scale="free_y") +
  theme(legend.position="bottom")

# ENSG00000154370.16 example,
test.UpwardRespondingGene.dat <- expression.dat %>%
  filter(Geneid == "ENSG00000154370.16") %>%
  mutate(treatment = factor(treatment))


fit1_up <- drm(formula = log2FC ~ doseInRisdiscale,
                 data = test.UpwardRespondingGene.dat,
                 fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,NA,0,NA)),
                 curveid = treatment,
                 pmodels=data.frame(treatment, 1, treatment, treatment),
                 robust = "mean"
)

fit2_up <- drm(formula = log2FC ~ doseInRisdiscale,
               data = test.UpwardRespondingGene.dat,
               fct = LL.3(names=c("Steepness", "UpperLimit", "ED50")),
               curveid = treatment,
               pmodels=data.frame(treatment, 1, treatment),
               robust = "mean"
)

fit3_up <- drm(formula = log2FC ~ doseInRisdiscale,
               data = test.UpwardRespondingGene.dat,
               fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,0,NA,NA)),
               curveid = treatment,
               pmodels=data.frame(treatment, treatment, 1, treatment),
               robust = "mean"
)

plot(fit1_up)
plot(fit2_up)
plot(fit3_up)

summary(fit1_up)
summary(fit2_up)
summary(fit3_up)

# ENSG00000115415.20 example
test.DownwardRespondingGene.dat <- expression.dat %>%
  filter(Geneid == "ENSG00000115415.20") %>%
  mutate(treatment = factor(treatment))

fit1_down <- drm(formula = log2FC ~ doseInRisdiscale,
                 data = test.DownwardRespondingGene.dat,
                 fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,NA,0,NA)),
                 curveid = treatment,
                 pmodels=data.frame(treatment, 1, treatment, treatment),
                 robust = "mean"
)

fit2_down <- drm(formula = log2FC ~ doseInRisdiscale,
               data = test.DownwardRespondingGene.dat,
               fct = LL.3(names=c("Steepness", "UpperLimit", "ED50")),
               curveid = treatment,
               pmodels=data.frame(treatment, 1, treatment),
               robust = "mean"
)

fit3_down <- drm(formula = log2FC ~ doseInRisdiscale,
               data = test.DownwardRespondingGene.dat,
               fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50"), fixed=c(NA,0,NA,NA)),
               curveid = treatment,
               pmodels=data.frame(treatment, treatment, 1, treatment),
               robust = "mean"
)


plot(fit1_down)
plot(fit2_down)
plot(fit3_down)

summary(fit1_down)
summary(fit2_down)
summary(fit3_down)



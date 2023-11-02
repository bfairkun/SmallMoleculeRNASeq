library(tidyverse)

dat <- read_tsv("config/samples.3MoleculesOfInterestExperiment.tsv") %>%
  distinct(SampleName) %>%
  separate(SampleName, into=c("Exp", "treatment", "dose.nM", "rep"), sep="_", remove=F)

dat.controls <- dat %>%
  filter(treatment == "DMSO")

dat.treatments <- dat %>%
  filter(!treatment == "DMSO") %>%
  group_by(treatment, dose.nM) %>%
  mutate(group_name = paste0(treatment, "_", dose.nM)) %>%
  ungroup()


dat.controls.expanded <- dat.controls %>%
  dplyr::select(SampleName) %>%
  expand(group_name = dat.treatments$group_name, nesting(SampleName))

# dat.treatments %>%
#   distinct(group_name) %>%
#   write_tsv("config/samples.3MoleculesOfInterestExperiment.Contrasts.tsv")

bind_rows(dat.treatments, dat.controls.expanded) %>%
  arrange(group_name, treatment) %>%
  dplyr::select(SampleName, group_name) %>%
  separate(SampleName, into=c("Exp", "treatment", "dose.nM", "rep"), sep="_", remove=F) %>%
  mutate(TreatmentOrControl = if_else(treatment == "DMSO", "DMSO", "treatment")) %>%
  dplyr::select(SampleName, TreatmentOrControl, group_name) %>%
  mutate(TreatmentOrControl = relevel(factor(TreatmentOrControl), "DMSO")) %>%
  arrange(group_name, TreatmentOrControl) %>%
  group_by(group_name) %>%
  group_walk(~ write_tsv(.x, paste0("SplicingAnalysis/Exp202310_3Molecules_Contrasts/",.y$group_name, ".groups.tsv"), col_names=F))

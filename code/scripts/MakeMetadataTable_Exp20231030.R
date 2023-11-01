library(tidyverse)

d <- read_tsv("config/samples.titrationseries.tsv")

dat <- read_tsv("scratch/NewSampleList.tsv", col_names = c("R1", "R2"))

dat %>%
  mutate(SampleIndentifier = str_replace(R1,"^.+?Sample_(.+?)_L00[12].+$", "\\1")) %>%
  mutate(SampleName = case_when(
    str_detect(SampleIndentifier, "^DMSO") ~ str_replace(SampleIndentifier, "DMSO-(.+?)_090823","Exp090823 DMSO NA \\1"),
    TRUE ~ str_replace(SampleIndentifier, "^(.+?)_(.+?)uM_(.+?$)","Exp090823 \\1 \\2 \\3"),
  )) %>%
  separate(SampleName, into=c("Experiment", "Treatment", "dose.nM", "Rep"), sep=' ') %>%
  mutate(dose.nM = case_when(
    dose.nM == "0-6" ~ 600,
    dose.nM == "0-12" ~ 120,
    TRUE ~ as.numeric(dose.nM) * 1000
  )) %>%
  unite(SampleName, Experiment, Treatment, dose.nM, Rep, remove = F, sep = "_") %>%
  dplyr::select(SampleName, R1, R2) %>%
  mutate(R1Hash = str_replace(R1, "^.+?_ds\\.(.+?)/.+$", "\\1")) %>%
  mutate(R2Hash = str_replace(R2, "^.+?_ds\\.(.+?)/.+$", "\\1")) %>%
  # dplyr::select(R1Hash, R2Hash)
  filter(R1Hash == R2Hash) %>%
  dplyr::select(SampleName, R1, R2) %>%
  write_tsv("config/samples.3MoleculesOfInterestExperiment.tsv")

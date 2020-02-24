#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# Data input
# Three datasets: GOBO, SCAN-B, TCGA
# Information per dataset: patient annotation table (all), gene expression data (all), gene information table (GOBO, TCGA)
#
#

# LOAD PACKAGES  -------------------------------------------------

library(tidyverse)
library(yaml)
library(Biobase)


# LOAD FUNCTIONS  -------------------------------------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# (1) PATIENT ANNOTATION   --------------------------------------------
# NOTE: IF YOU ONLY WANT THE IDS IN A TABLE, GO STRAIGHT TO (1A)

# GOBO
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")
setwd(config$input_file_paths$directory)
patient_annotation_gobo <- loadRData(config$input_file_paths$patient_annotation)

# SCAN-B
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_scanb.yml")
setwd(config$input_file_paths$directory)
patient_annotation_scanb <- loadRData(config$input_file_paths$patient_annotation)

# TCGA
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")
setwd(config$input_file_paths$directory)
patient_annotation_tcga <- pData(readRDS(config$input_file_paths$rds))


# to create a table with only ids for all datasets
# GOBO
patient_annotation_gobo <- patient_annotation_gobo["SampleID"]
patient_annotation_gobo <-patient_annotation_gobo %>% rename(sample.id = SampleID)
patient_annotation_gobo$dataset <- "GOBO"

# SCAN-B
patient_annotation_scanb <- patient_annotation_scanb["rba"]
patient_annotation_scanb <- patient_annotation_scanb %>% rename(sample.id = rba)
patient_annotation_scanb$dataset <- "SCAN-B"

# TCGA
patient_annotation_tcga <- patient_annotation_tcga["sample_id"]
patient_annotation_tcga <- patient_annotation_tcga %>% rename(sample.id = sample_id)
patient_annotation_tcga$dataset <- "TCGA"

patient_ids_all_datasets <- bind_rows(patient_annotation_tcga, patient_annotation_gobo, patient_annotation_scanb)

# write.csv(patient_ids_all_datasets, "/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_only_ids.csv", row.names=FALSE)

rm(patient_annotation_tcga)
rm(patient_annotation_gobo)
rm(patient_annotation_scanb)


# (1A) ONLY PATIENT IDS FROM ANNOTATION   -----------------------------------

patient_ids_all_datasets <- read.csv("/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_only_ids.csv")



# (2) GENE TABLE INFORMATION  ---------------------------------

# GOBO
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")
setwd(config$input_file_paths$directory)
gene_table_gobo <- loadRData(config$input_file_paths$gene_table)

# SCAN-B
# none

# TCGA
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")
setwd(config$input_file_paths$directory)
gene_table_tcga <- fData(readRDS(config$input_file_paths$rds))

rm(config)



# (3) GEX  ----------------------------------------------------

# GOBO
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")
setwd(config$input_file_paths$directory)
gex_matrix_gobo <- loadRData(config$input_file_paths$gex_matrix)

# SCAN-B
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_scanb.yml")
setwd(config$input_file_paths$directory)
gex_matrix_scanb <- loadRData(config$input_file_paths$gex_matrix)

# TCGA
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")
setwd(config$input_file_paths$directory)
gex_matrix_tcga <- exprs(readRDS(config$input_file_paths$rds))

rm(config)



# CLAMS RESULTS  ----------------------------------------------------

patient_annotation_clams <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/0_clams/clams_all_samples.csv")



# IMMUNE ANALYSIS RESULTS  -------------------------------------------

# list of 65 genes from immune signature module that are present in all 3 datasets (GOBO, SCAN-B, TCGA)
immune_genes_in_common <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_genes_in_common.txt")
immune_genes_in_common <- immune_genes_in_common$V1

# table of patient ids and immune signature values (sum of ranks) for all 3 datasets
# sum = sum FPKM when gene symbol appears more than once, max = get the highest value when gene symbol is there more than once
patient_annotation_immune <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_samples.csv")

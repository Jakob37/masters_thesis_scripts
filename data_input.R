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


# (1A) to create a table with ids, dataset and cancer type for all datasets
# GOBO
patient_annotation_gobo <- patient_annotation_gobo["SampleID"]
patient_annotation_gobo <-patient_annotation_gobo %>% rename(sample.id = SampleID)
patient_annotation_gobo$cancer.type <- "BRCA"
patient_annotation_gobo$dataset <- "GOBO"

# SCAN-B
patient_annotation_scanb <- patient_annotation_scanb["rba"]
patient_annotation_scanb <- patient_annotation_scanb %>% rename(sample.id = rba)
patient_annotation_scanb$cancer.type <- "BRCA"
patient_annotation_scanb$dataset <- "SCAN-B"

# TCGA
patient_annotation_tcga <- patient_annotation_tcga[c("sample_id", "cancer.type")]
patient_annotation_tcga <- patient_annotation_tcga %>% rename(sample.id = sample_id)
patient_annotation_tcga$dataset <- "TCGA"

patient_ids_all_datasets <- bind_rows(patient_annotation_tcga, patient_annotation_gobo, patient_annotation_scanb)

# write.csv(patient_ids_all_datasets, "/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_only_ids.csv", row.names=FALSE)


# (1B) to create a table with os, dataset and cancer type for all datasets
# GOBO
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")
patient_annotation_gobo <- patient_annotation_gobo[c(config$patient_annotation_info$sample_column_name,
                                                     config$patient_annotation_info$overall_survival_time_column_name,
                                                     config$patient_annotation_info$overall_survival_event_column_name)]
patient_annotation_gobo <-patient_annotation_gobo %>% rename(sample.id = config$patient_annotation_info$sample_column_name,
                                                             OS.time = config$patient_annotation_info$overall_survival_time_column_name,
                                                             OS.event = config$patient_annotation_info$overall_survival_event_column_name)
patient_annotation_gobo$OS.time.type <- config$patient_annotation_info$overall_survival_time_column_unit_table
patient_annotation_gobo$cancer.type <- "BRCA"
patient_annotation_gobo$dataset <- "GOBO"

# SCAN-B
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_scanb.yml")
patient_annotation_scanb <- patient_annotation_scanb[c(config$patient_annotation_info$sample_column_name,
                                                     config$patient_annotation_info$overall_survival_time_column_name,
                                                     config$patient_annotation_info$overall_survival_event_column_name)]
patient_annotation_scanb <-patient_annotation_scanb %>% rename(sample.id = config$patient_annotation_info$sample_column_name,
                                                             OS.time = config$patient_annotation_info$overall_survival_time_column_name,
                                                             OS.event = config$patient_annotation_info$overall_survival_event_column_name)
patient_annotation_scanb$OS.time.type <- config$patient_annotation_info$overall_survival_time_column_unit_table
patient_annotation_scanb$cancer.type <- "BRCA"
patient_annotation_scanb$dataset <- "SCAN-B"

# TCGA
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")
patient_annotation_tcga <- patient_annotation_tcga[c(config$patient_annotation_info$sample_column_name,
                                                       config$patient_annotation_info$overall_survival_time_column_name,
                                                       config$patient_annotation_info$overall_survival_event_column_name,
                                                       config$patient_annotation_info$cancer_type_column)]
patient_annotation_tcga <-patient_annotation_tcga %>% rename(sample.id = config$patient_annotation_info$sample_column_name,
                                                               OS.time = config$patient_annotation_info$overall_survival_time_column_name,
                                                               OS.event = config$patient_annotation_info$overall_survival_event_column_name)
patient_annotation_tcga$OS.time.type <- config$patient_annotation_info$overall_survival_time_column_unit_table
patient_annotation_tcga$dataset <- "TCGA"

patient_ids_all_datasets <- bind_rows(patient_annotation_tcga, patient_annotation_gobo, patient_annotation_scanb)
factor_cols <- c("OS.event", "cancer.type", "OS.time.type", "dataset")
patient_ids_all_datasets[factor_cols] <- lapply(patient_ids_all_datasets[factor_cols], factor)

# write.csv(patient_ids_all_datasets, "/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_os.csv", row.names=FALSE)

rm(patient_annotation_tcga)
rm(patient_annotation_gobo)
rm(patient_annotation_scanb)


# (1A) ONLY PATIENT IDS FROM ANNOTATION   -----------------------------------

patient_ids_all_datasets <- read.csv("/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_only_ids.csv")
patient_ids_all_datasets$sample.id <- as.character(patient_ids_all_datasets$sample.id)


# (1B) OS FROM ANNOTATION   -----------------------------------

patient_os_all_datasets <- read.csv("/media/deboraholi/Data/LUND/9 THESIS/data/all_samples_os.csv")


# (2) GENE TABLE INFORMATION  ---------------------------------

# GOBO
setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")
setwd(config$input_file_paths$directory)
gene_table_gobo <- loadRData(config$input_file_paths$gene_table)

# SCAN-B
gene_table_scanb <- loadRData("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/JVCs.predictors_Helena/Gene.ID.Entrez.Rdata") # from Helena (JVC.predictor)

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


# PLOTS NAMES  ------------------------------------------------------

source("/media/deboraholi/Data/LUND/9 THESIS/src/plots_names.R")



# CLAMS RESULTS  ----------------------------------------------------

patient_annotation_clams <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/0_clams/clams_all_samples.csv")



# PROLIFERATION ANALYSIS RESULTS  -----------------------------------

# list of BLAH genes from proliferation module that are present in all 3 datasets (GOBO, SCAN-B, TCGA)
prolif_fred_genes_in_common <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/prolif_fred_genes_in_common.txt")
prolif_fred_genes_in_common <- prolif_fred_genes_in_common$V1

# list of BLAH genes from proliferation module that are present in all 3 datasets (GOBO, SCAN-B, TCGA)
prolif_karl_genes_in_common <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/prolif_karl_genes_in_common.txt")
prolif_karl_genes_in_common <- prolif_karl_genes_in_common$V1

# table of patient ids, clams.class, proliferation values (sum of ranks) for all 3 datasets (fred+karl), karl groups
patient_annotation_prolif <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/prolif_all_samples.csv")



# IMMUNE ANALYSIS RESULTS  -------------------------------------------

# list of 65 genes from immune signature module that are present in all 3 datasets (GOBO, SCAN-B, TCGA)
immune_genes_in_common <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_genes_in_common.txt")
immune_genes_in_common <- immune_genes_in_common$V1

# table of patient ids and immune signature values (sum of ranks) for all 3 datasets
# sum = sum FPKM when gene symbol appears more than once, max = get the highest value when gene symbol is there more than once
patient_annotation_immune <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_samples.csv")



# FUNCTION TO RUN DIFFERENT SSPS  --------------------------------------

source("/media/deboraholi/Data/LUND/9 THESIS/src/run_ssps_function.R")



# ROR SSPs  -------------------------------------------------------------

# reduced version
ror.red.aims.gs <- loadRData("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/Training_Run9564Genes_noNorm_SSP.scaled.ROR.tot.asT0.c005.Fcc15_5x5foldCV.num.rules.50_50.selRules.AIMS.GS.RData")

# complete version
ror.all.aims.gs <- loadRData("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/Training_Run19081Genes_noNorm_SSP.scaled.ROR.tot.asT0.c005.Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData")



# ROR RESULTS  ------------------------------------------------------------

patient_annotation_ror <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ror_all_samples.csv")



# PAM50 SSPs  -------------------------------------------------------------

# reduced version
load("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/Training_Run9564Genes_noNorm_SSP.PAM50subtype4Most.Fcc15_5x5foldCV.num.rules.50_22.selRules.AIMS.GS.RData")

# complete version
load("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/Training_Run19081Genes_noNorm_SSP.PAM50subtype4Most.Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData")



# PAM50 RESULTS  ------------------------------------------------------------

patient_annotation_pam50 <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/pam50_all_samples.csv")
#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO PREPARE BREAST CANCER DATA, RUN ROR SSP, ANALYSE DATA
#
#


# CLEAR ENVIRONMENT ----------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX
# FUNCTION TO RUN DIFFERENT SSPS
# ROR SSPs

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/")


# ANALYSIS  ------------------------------------------------------

# Functions
get_sample_class_from_predictor_result <- function(ror_result) {
  sample.class <- ror_result$cl[,]
  sample.class <- data.frame(as.list(sample.class))
  sample.class <- t(sample.class)
  sample.class <- rownames_to_column(data.frame(sample.class), var = "sample.id")
}


# GOBO

# Create an entrez id list from the reporter ids present in the gex_matrix  
entrez_ids_gobo <- gex_matrix_gobo[,1:2]
entrez_ids_gobo <- rownames_to_column(data.frame(entrez_ids_gobo), var = "reporterId")
entrez_ids_gobo <- left_join(entrez_ids_gobo, gene_table_gobo[,c("reporterId","entrezId")])
entrez_ids_gobo <- entrez_ids_gobo %>% mutate(entrezId = paste0("e",entrezId))
entrez_ids_gobo <- entrez_ids_gobo$entrezId

rm(gene_table_gobo)

# Rename gex_matrix rows with entrez ids instead of reporter ids
row.names(gex_matrix_gobo) <- entrez_ids_gobo

# run ROR predictor
ror_red_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, ror.red.aims.gs)
ror_all_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_gobo <- get_sample_class_from_predictor_result(ror_red_result_gobo)
ror_red_class_gobo <- ror_red_class_gobo %>% rename(ror.red.class = sample.class)
ror_all_class_gobo <- get_sample_class_from_predictor_result(ror_all_result_gobo)
ror_all_class_gobo <- ror_all_class_gobo %>% rename(ror.all.class = sample.class)

# save to patient id table
ror_class_gobo <- inner_join(patient_ids_all_datasets, ror_red_class_gobo, by="sample.id")
ror_class_gobo <- inner_join(ror_class_gobo, ror_all_class_gobo, by="sample.id")

rm(gex_matrix_gobo)
rm(entrez_ids_gobo)


# SCAN-B

# Create an entrez id list from the gene symbols present in the gex_matrix  
entrez_ids_scanb <- gex_matrix_scanb[,1:2]
entrez_ids_scanb <- rownames_to_column(data.frame(entrez_ids_scanb), var = "geneSymbol")
entrez_ids_scanb <- left_join(entrez_ids_scanb, gene_table_scanb[,c("Gene.Name","EntrezGene")], by=c("geneSymbol" = "Gene.Name"))
entrez_ids_scanb <- entrez_ids_scanb %>% mutate(entrezId = paste0("e",EntrezGene))
entrez_ids_scanb <- entrez_ids_scanb$entrezId

rm(gene_table_gobo)

# Rename gex_matrix rows with entrez ids instead of gene symbols
row.names(gex_matrix_scanb) <- entrez_ids_scanb

# run ROR predictor
ror_red_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, ror.red.aims.gs)
ror_all_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_scanb <- get_sample_class_from_predictor_result(ror_red_result_scanb)
ror_red_class_scanb <- ror_red_class_scanb %>% rename(ror.red.class = sample.class)
ror_all_class_scanb <- get_sample_class_from_predictor_result(ror_all_result_scanb)
ror_all_class_scanb <- ror_all_class_scanb %>% rename(ror.all.class = sample.class)

# save to patient id table
ror_class_scanb <- inner_join(patient_ids_all_datasets, ror_red_class_scanb, by="sample.id")
ror_class_scanb <- inner_join(ror_class_scanb, ror_all_class_scanb, by="sample.id")

rm(gex_matrix_scanb)
rm(entrez_ids_scanb)


# TCGA

# Create an entrez id list from the gene symbols present in the gex_matrix  
entrez_ids_tcga <- gex_matrix_tcga[,1:2]
entrez_ids_tcga <- rownames_to_column(data.frame(entrez_ids_tcga), var = "reporterId")
entrez_ids_tcga <- left_join(entrez_ids_tcga, gene_table_tcga[,c("ENSG","ENTREZID")], by=c("reporterId" = "ENSG"))
entrez_ids_tcga <- entrez_ids_tcga %>% mutate(entrezId = paste0("e",ENTREZID))
entrez_ids_tcga <- entrez_ids_tcga$entrezId

rm(gene_table_tcga)

# Rename gex_matrix rows with entrez ids instead of reporter ids
row.names(gex_matrix_tcga) <- entrez_ids_tcga

# run ROR predictor
cancer_types <- levels(patient_ids_all_datasets$cancer.type)
for (cancer_type in cancer_types) {
  
  # Subset data according to cancer type
  print(cancer_type)
  sample_list <- subset(patient_ids_all_datasets, dataset == 'TCGA' & cancer.type == cancer_type) %>% pull(sample.id)
  col_num <- which(colnames(gex_matrix_tcga) %in% sample_list)
  new_gex_matrix <- gex_matrix_tcga[,c(col_num)]
  gene_entrez <- row.names(new_gex_matrix)
  
  # run ROR predictors
  print('Running reduced')
  ror_red_result_tcga <- applyAIMS(new_gex_matrix, gene_entrez, ror.red.aims.gs)
  print('Running complete')
  ror_all_result_tcga <- applyAIMS(new_gex_matrix, gene_entrez, ror.all.aims.gs)
  
  # extract the classification result
  if (exists("ror_red_class_tcga")) {
    ror_red_class_tcga <- bind_rows(ror_red_class_tcga, get_sample_class_from_predictor_result(ror_red_result_tcga)) 
    } else {
      ror_red_class_tcga <- get_sample_class_from_predictor_result(ror_red_result_tcga)
    }
  
  if (exists("ror_all_class_tcga")) {
    ror_all_class_tcga <- bind_rows(ror_all_class_tcga, get_sample_class_from_predictor_result(ror_all_result_tcga))
    } else {
      ror_all_class_tcga <- get_sample_class_from_predictor_result(ror_all_result_tcga)
    }
  
  rm(new_gex_matrix)
}

ror_red_class_tcga <- ror_red_class_tcga %>% rename(ror.red.class = sample.class)
ror_all_class_tcga <- ror_all_class_tcga %>% rename(ror.all.class = sample.class)

# save to patient id table
ror_class_tcga <- inner_join(patient_ids_all_datasets, ror_red_class_tcga, by="sample.id")
ror_class_tcga <- inner_join(ror_class_tcga, ror_all_class_tcga, by="sample.id")

rm(gex_matrix_tcga)
rm(entrez_ids_tcga)

# combine all datasets
patient_annotation_ror <- bind_rows(ror_class_gobo, ror_class_scanb, ror_class_tcga)
patient_annotation_ror$ror.red.class <- factor(patient_annotation_ror$ror.red.class)
patient_annotation_ror$ror.all.class <- factor(patient_annotation_ror$ror.all.class)

# RESULT
# Save final table
write.csv(patient_annotation_ror, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ror_all_samples.csv", row.names=FALSE)



# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# ROR RESULTS 
# CLAMS RESULTS
# PLOTS NAMES

# compare classification between predictors
table(patient_annotation_ror$ror.red.class, patient_annotation_ror$ror.all.class)

# compare what is low (<40) and what is High (>60)

# compare with CLAMS class

# do OS?


#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO PREPARE CANCER DATA, RUN CLAMS PREDICTOR, ANALYSE DATA
#
#


# CLEAR ENVIRONMENT ----------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(yaml)



# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX



# CLAMS ANALYSIS  -----------------------------------

library(CLAMS)
library(e1071)


# GOBO

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")

# Create a gene symbol list from the reporter ids present in the gex_matrix  
gene_symbols_gobo <- gex_matrix_gobo[,1:2]
gene_symbols_gobo <- rownames_to_column(data.frame(gene_symbols_gobo), var = "reporterId")
gene_symbols_gobo <- left_join(gene_symbols_gobo, gene_table_gobo[,c("reporterId","geneSymbol")])
gene_symbols_gobo <- gene_symbols_gobo$geneSymbol

rm(gene_table_gobo)

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_gobo) <- gene_symbols_gobo

# run CLAMS predictor
clams_result_gobo <- applyCLAMS(gex_matrix_gobo, gene_symbols_gobo)
# extract the classification result from CLAMS
clams.class <- clams_result_gobo$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_gobo <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_gobo)
rm(gene_symbols_gobo)



# SCAN-B

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_scanb.yml")

# Create a gene symbol list from gex_matrix
gene_symbols_scanb <- gex_matrix_scanb[,1:2]
gene_symbols_scanb <- rownames_to_column(data.frame(gex_matrix_scanb), var = "geneSymbol")
gene_symbols_scanb <- gene_symbols_scanb$geneSymbol

# run CLAMS predictor
clams_result_scanb <- applyCLAMS(gex_matrix_scanb, gene_symbols_scanb)
# extract the classification result from CLAMS
clams.class <- clams_result_scanb$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_scanb <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_scanb)
rm(gene_symbols_scanb)


# TCGA

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")

# Create a gene symbol list from the reporter ids present in the gex_matrix  
gene_symbols_tcga <- gex_matrix_tcga[,1:2]
gene_symbols_tcga <- rownames_to_column(data.frame(gene_symbols_tcga), var = "reporterId")
gene_symbols_tcga <- left_join(gene_symbols_tcga, gene_table_tcga[,c("ENSG", "SYMBOL")], by=c("reporterId" = "ENSG"))
gene_symbols_tcga <- gene_symbols_tcga$SYMBOL

rm(gene_table_tcga)

# Replace VEGFD for FIGF
gene_symbols_tcga <- as.list(gene_symbols_tcga)
gene_symbols_tcga <- rapply(gene_symbols_tcga, function(x) ifelse(x=="VEGFD","FIGF",x), how="replace")

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_tcga) <- gene_symbols_tcga

rm(gene_symbols_tcga)

# run CLAMS predictor for subsets of each cancer type
cancer_types <- levels(factor(patient_annotation_tcga$cancer.type))
for (cancer_type in cancer_types) {
  
  # Subset data according to cancer type
  print(cancer_type)
  sample_list <- subset(patient_annotation_tcga, cancer.type == cancer_type) %>% pull(sample_id)
  col_num <- which(colnames(gex_matrix_tcga) %in% sample_list)
  new_gex_matrix <- gex_matrix_tcga[,c(col_num)]
  gene_symbols <- row.names(new_gex_matrix)
  
  # run CLAMS predictor
  clams_result <- applyCLAMS(new_gex_matrix, gene_symbols)
  # extract the classification result from CLAMS
  clams.class <- clams_result$cl[,]
  
  clams.class <- data.frame(as.list(clams.class))
  clams.class <- t(clams.class)
  clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")
  
  if (exists("clams_class_tcga")) {
    clams_class_tcga <- bind_rows(clams_class_tcga, clams.class) } else {
      clams_class_tcga <- clams.class
    }
  
  rm(new_gex_matrix)
}

clams_class_tcga <- inner_join(patient_ids_all_datasets, clams_class_tcga, by="sample.id")

rm(patient_annotation_tcga)
rm(gex_matrix_tcga)


# combine all datasets
patient_annotation_clams <- bind_rows(clams_class_gobo, clams_class_scanb, clams_class_tcga)
# transmute bronchioid to TRU and nonbronchioid to NonTRU
patient_annotation_clams <- patient_annotation_clams %>% 
                              mutate(clams.class = case_when(clams.class == "bronchioid" ~ "TRU",
                                     clams.class == "nonbronchioid" ~ "NonTRU",
                                     TRUE ~ clams.class))
patient_annotation_clams$clams.class <- factor(patient_annotation_clams$clams.class, levels=c("TRU", "NonTRU"))

  
# RESULT
# Save final table
write.csv(patient_annotation_clams, "/media/deboraholi/Data/LUND/9 THESIS/0_clams/clams_all_samples.csv", row.names=FALSE)


  
# CHANGE FROM HERE ------------------



# GENERAL VISUALIZATION OF PATIENT ANNOTATION INFORMATION BY CLAMS CLASS (ONE TYPE IN DATASET) --------------------------------------------------

# Barplots
if (length(config$see_variables_by_clams$barplot) > 0) {
  for (i in 1:length(config$see_variables_by_clams$barplot)) {
    variable <- config$see_variables_by_clams$barplot[i]
    this_plot <- ggplot(data=subset(patient_annotation, !is.na(get(variable))), aes_string(x=variable)) +
      geom_bar(aes(fill=clams_class)) +
      labs(y ="Number of samples") +
      theme_minimal() +
      scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}


# Histograms
if (length(config$see_variables_by_clams$histogram) > 0) {
  for (i in 1:length(config$see_variables_by_clams$histogram)) {
    variable <- config$see_variables_by_clams$histogram[i]
    this_plot <- ggplot(data=patient_annotation, aes_string(x=variable)) +
      geom_histogram(aes(fill=clams_class)) +
      labs(y ="Number of samples") +
      theme_minimal() +
      scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}


# Frequency polygons
if (length(config$see_variables_by_clams$histogram) > 0) {
  for (i in 1:length(config$see_variables_by_clams$histogram)) {
    variable <- config$see_variables_by_clams$histogram[i]
    this_plot <- ggplot(data=patient_annotation, aes(get(variable), stat(density), colour=clams_class)) +
      geom_freqpoly(bins=15) +
      labs(x = variable) +
      theme_minimal() +
      scale_colour_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}

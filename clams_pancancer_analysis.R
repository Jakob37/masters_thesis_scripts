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
library(org.Hs.eg.db)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX
# functions.R


# CLAMS ANALYSIS  -----------------------------------

library(CLAMS)
library(e1071)


# GOBO

# Create a gene symbol list from the reporter ids present in the gex_matrix  
gene_symbols_gobo <- gex_matrix_gobo[,1:2]
gene_symbols_gobo <- rownames_to_column(data.frame(gene_symbols_gobo), var = "reporterId")
gene_symbols_gobo <- left_join(gene_symbols_gobo, gene_table_gobo[,c("reporterId","geneSymbol")])
gene_symbols_gobo <- gene_symbols_gobo$geneSymbol

rm(gene_table_gobo)

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_gobo) <- gene_symbols_gobo

# get CLAMS genes and filter gex_matrix for only those genes
clams_genes <- get.all.pairs.genes(CLAMSmodel$all.pairs)
gex_matrix_gobo_clams <- gex_matrix_gobo[clams_genes,]

# run CLAMS predictor
clams_result_gobo <- applyCLAMS(gex_matrix_gobo_clams, row.names(gex_matrix_gobo_clams))

# extract the classification result from CLAMS
clams.class <- clams_result_gobo$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_gobo <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_gobo)
rm(gene_symbols_gobo)


# SCAN-B

# Create a gene symbol list from gex_matrix
gene_symbols_scanb <- row.names(gex_matrix_scanb)

# get CLAMS genes and filter gex_matrix for only those genes
clams_genes <- get.all.pairs.genes(CLAMSmodel$all.pairs)
gex_matrix_scanb_clams <- gex_matrix_scanb[clams_genes,]

# run CLAMS predictor
clams_result_scanb <- applyCLAMS(gex_matrix_scanb_clams, row.names(gex_matrix_scanb_clams))

# extract the classification result from CLAMS
clams.class <- clams_result_scanb$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_scanb <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_scanb)
rm(gene_symbols_scanb)


# TCGA

# Create an entrez id list from the ensembl ids present in the gex_matrix
ensembl_ids_tcga <- row.names(gex_matrix_tcga)
gene_symbols_tcga <- mapIds(org.Hs.eg.db, ensembl_ids_tcga, 'SYMBOL', 'ENSEMBL')
gene_symbols_tcga <- unname(gene_symbols_tcga)

# Replace VEGFD for FIGF
gene_symbols_tcga <- as.list(gene_symbols_tcga)
gene_symbols_tcga <- rapply(gene_symbols_tcga, function(x) ifelse(x=="VEGFD","FIGF",x), how="replace")

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_tcga) <- gene_symbols_tcga

# get CLAMS genes and filter gex_matrix for only those genes
clams_genes <- get.all.pairs.genes(CLAMSmodel$all.pairs)
gex_matrix_tcga_clams <- gex_matrix_tcga[clams_genes,]

# run CLAMS predictor
clams_result_tcga <- applyCLAMS(gex_matrix_tcga_clams, row.names(gex_matrix_tcga_clams))

# extract the classification result from CLAMS
clams.class <- clams_result_tcga$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_tcga <- inner_join(patient_ids_all_datasets, clams_class_tcga, by="sample.id")

# for only BRCA
brca_tcga_subset_samples <- subset(patient_ids_all_datasets, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
brca_tcga_gex_matrix <- gex_matrix_tcga_filtered[,brca_tcga_subset_samples]

clams_result_brca_tcga <- applyCLAMS(brca_tcga_gex_matrix, row.names(brca_tcga_gex_matrix))


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



## LIHC  ---------------------------------------------------

# load lihc_subtypes_huang_2020

lihc <- subset(patient_annotation_tcga, cancer.type == 'LIHC')
lihc_subtypes_huang_2020 <- left_join(lihc, lihc_subtypes_huang_2020, by='sample_barcode')
lihc_subtypes_huang_2020 <- left_join(lihc_subtypes_huang_2020, patient_annotation_clams[c("sample.id", "clams.class")], 
                                      by=c('sample_id'='sample.id'))

addmargins(table(lihc_subtypes_huang_2020$clams.class, lihc_subtypes_huang_2020$subtype, useNA='ifany'))
#          1   2   3   4 <NA> Sum
# NonTRU  77 101  95  55    2 330
# TRU      9   1  25   3    0  38
# Sum     86 102 120  58    2 368


# NEEDS TO BE CHANGED FROM HERE ------  
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

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
# PAM50 SSPs

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/")


# ANALYSIS  ------------------------------------------------------

# Functions
get_sample_class_from_predictor_result <- function(pam50_result) {
  sample.class <- pam50_result$cl[,]
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

# run PAM50 predictor
pam50_red_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, pam50.red.aims.gs)
pam50_all_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, pam50.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
pam50_red_class_gobo <- get_sample_class_from_predictor_result(pam50_red_result_gobo)
pam50_red_class_gobo <- pam50_red_class_gobo %>% rename(pam50.red.class = sample.class)
pam50_all_class_gobo <- get_sample_class_from_predictor_result(pam50_all_result_gobo)
pam50_all_class_gobo <- pam50_all_class_gobo %>% rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_gobo <- inner_join(patient_ids_all_datasets, pam50_red_class_gobo, by="sample.id")
pam50_class_gobo <- inner_join(pam50_class_gobo, pam50_all_class_gobo, by="sample.id")

rm(gex_matrix_gobo)
rm(entrez_ids_gobo)


# SCAN-B

# Create an entrez id list from the gene symbols present in the gex_matrix  
entrez_ids_scanb <- gex_matrix_scanb[,1:2]
entrez_ids_scanb <- rownames_to_column(data.frame(entrez_ids_scanb), var = "geneSymbol")
entrez_ids_scanb <- left_join(entrez_ids_scanb, gene_table_scanb[,c("Gene.Name","EntrezGene")], by=c("geneSymbol" = "Gene.Name"))
entrez_ids_scanb <- entrez_ids_scanb %>% mutate(entrezId = paste0("e",EntrezGene))
entrez_ids_scanb <- entrez_ids_scanb$entrezId

rm(gene_table_scanb)

# Rename gex_matrix rows with entrez ids instead of gene symbols
row.names(gex_matrix_scanb) <- entrez_ids_scanb

# run PAM50 predictor
pam50_red_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, pam50.red.aims.gs)
pam50_all_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, pam50.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
pam50_red_class_scanb <- get_sample_class_from_predictor_result(pam50_red_result_scanb)
pam50_red_class_scanb <- pam50_red_class_scanb %>% rename(pam50.red.class = sample.class)
pam50_all_class_scanb <- get_sample_class_from_predictor_result(pam50_all_result_scanb)
pam50_all_class_scanb <- pam50_all_class_scanb %>% rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_scanb <- inner_join(patient_ids_all_datasets, pam50_red_class_scanb, by="sample.id")
pam50_class_scanb <- inner_join(pam50_class_scanb, pam50_all_class_scanb, by="sample.id")

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

# run PAM50 predictor
cancer_types <- levels(patient_ids_all_datasets$cancer.type)
for (cancer_type in cancer_types) {
  
  # Subset data according to cancer type
  print(cancer_type)
  sample_list <- subset(patient_ids_all_datasets, dataset == 'TCGA' & cancer.type == cancer_type) %>% pull(sample.id)
  col_num <- which(colnames(gex_matrix_tcga) %in% sample_list)
  new_gex_matrix <- gex_matrix_tcga[,c(col_num)]
  gene_entrez <- row.names(new_gex_matrix)
  
  # run PAM50 predictors
  print('Running reduced')
  pam50_red_result_tcga <- applyAIMS(new_gex_matrix, gene_entrez, pam50.red.aims.gs)
  print('Running complete')
  pam50_all_result_tcga <- applyAIMS(new_gex_matrix, gene_entrez, pam50.all.aims.gs)
  
  # extract the classification result
  if (exists("pam50_red_class_tcga")) {
    pam50_red_class_tcga <- bind_rows(pam50_red_class_tcga, get_sample_class_from_predictor_result(pam50_red_result_tcga)) 
    } else {
      pam50_red_class_tcga <- get_sample_class_from_predictor_result(pam50_red_result_tcga)
    }
  
  if (exists("pam50_all_class_tcga")) {
    pam50_all_class_tcga <- bind_rows(pam50_all_class_tcga, get_sample_class_from_predictor_result(pam50_all_result_tcga))
    } else {
      pam50_all_class_tcga <- get_sample_class_from_predictor_result(pam50_all_result_tcga)
    }
  
  rm(new_gex_matrix)
}

pam50_red_class_tcga <- pam50_red_class_tcga %>% rename(pam50.red.class = sample.class)
pam50_all_class_tcga <- pam50_all_class_tcga %>% rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_tcga <- inner_join(patient_ids_all_datasets, pam50_red_class_tcga, by="sample.id")
pam50_class_tcga <- inner_join(pam50_class_tcga, pam50_all_class_tcga, by="sample.id")

rm(gex_matrix_tcga)
rm(entrez_ids_tcga)

# combine all datasets
patient_annotation_pam50 <- bind_rows(pam50_class_gobo, pam50_class_scanb, pam50_class_tcga)
patient_annotation_pam50$pam50.red.class <- factor(patient_annotation_pam50$pam50.red.class)
patient_annotation_pam50$pam50.all.class <- factor(patient_annotation_pam50$pam50.all.class)

# RESULT
# Save final table
write.csv(patient_annotation_pam50, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/pam50_all_samples.csv", row.names=FALSE)



# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# PAM50 RESULTS 
# CLAMS RESULTS


# compare classification between predictors
table(patient_annotation_pam50$pam50.red.class, 
        patient_annotation_pam50$pam50.all.class, useNA='ifany')
prop.table(table(patient_annotation_pam50$pam50.red.class, 
                 patient_annotation_pam50$pam50.all.class, useNA='ifany'))


# with CLAMS classification as well
patient_annotation_pam50_clams <- left_join(patient_annotation_pam50, patient_annotation_clams, 
                                          by=c("sample.id", "cancer.type", "dataset"))
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.red.class, useNA='ifany')
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.all.class, useNA='ifany')

patient_annotation_pam50_clams %>%
  ggplot(aes(x=pam50.red.class, fill=clams.class)) +
  geom_histogram(position = 'identity', stat="count") +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1")) +
  theme_minimal() +
  facet_wrap(~clams.class) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  labs(y = "Number of samples", x = "PAM50 reduced classification")

ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/PAM50_reduced_clams.png", width=9.3, height=5, dpi=300)

patient_annotation_pam50_clams %>%
  ggplot(aes(x=pam50.all.class, fill=clams.class)) +
  geom_histogram(position = 'identity', stat="count") +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1")) +
  theme_minimal() +
  facet_wrap(~clams.class) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  labs(y = "Number of samples", x = "PAM50 not reduced classification")

ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/PAM50_not_reduced_clams.png", width=9.3, height=5, dpi=300)

# do OS?


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
library(viridisLite)

# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX
# FUNCTION TO RUN DIFFERENT SSPS run_ssps_function.R
# GENERAL FUNCTIONS: GET SAMPLES CLASS FROM AIMS RESULT functions.R
# ROR SSPs

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/")


# ANALYSIS  ------------------------------------------------------

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
ror_red_class_gobo <- get_samples_class_from_predictor_result(ror_red_result_gobo)
ror_red_class_gobo <- ror_red_class_gobo %>% rename(ror.red.class = sample.class)
ror_all_class_gobo <- get_samples_class_from_predictor_result(ror_all_result_gobo)
ror_all_class_gobo <- ror_all_class_gobo %>% rename(ror.all.class = sample.class)

# save to patient id table
ror_class_gobo <- inner_join(patient_ids_all_datasets, ror_red_class_gobo, by="sample.id")
ror_class_gobo <- left_join(ror_class_gobo, ror_all_class_gobo, by="sample.id")

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

# run ROR predictor
ror_red_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, ror.red.aims.gs)
ror_all_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_scanb <- get_samples_class_from_predictor_result(ror_red_result_scanb)
ror_red_class_scanb <- ror_red_class_scanb %>% rename(ror.red.class = sample.class)
ror_all_class_scanb <- get_samples_class_from_predictor_result(ror_all_result_scanb)
ror_all_class_scanb <- ror_all_class_scanb %>% rename(ror.all.class = sample.class)

# save to patient id table
ror_class_scanb <- inner_join(patient_ids_all_datasets, ror_red_class_scanb, by="sample.id")
ror_class_scanb <- left_join(ror_class_scanb, ror_all_class_scanb, by="sample.id")

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
  assign(paste0("ror_red_result_tcga_", cancer_type), ror_red_result_tcga)
  print('Running complete')
  ror_all_result_tcga <- applyAIMS(new_gex_matrix, gene_entrez, ror.all.aims.gs)
  assign(paste0("ror_all_result_tcga_", cancer_type), ror_all_result_tcga)
  
  #extract the classification result
  if (exists("ror_red_class_tcga")) {
    ror_red_class_tcga <- bind_rows(ror_red_class_tcga, get_samples_class_from_predictor_result(ror_red_result_tcga))
    } else {
      ror_red_class_tcga <- get_samples_class_from_predictor_result(ror_red_result_tcga)
    }
  if (exists("ror_all_class_tcga")) {
    ror_all_class_tcga <- bind_rows(ror_all_class_tcga, get_samples_class_from_predictor_result(ror_all_result_tcga))
    } else {
      ror_all_class_tcga  <- get_samples_class_from_predictor_result(ror_all_result_tcga)
    }
  rm(new_gex_matrix)
  rm(ror_red_result_tcga)
  rm(ror_all_result_tcga)
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
# write.csv(patient_annotation_ror, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ror_all_samples.csv", row.names=FALSE)



# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R
# ROR RESULTS 
# CLAMS RESULTS


# Pancancer  -------
# compare what is low (<40) and what is High (>60)
ror.low <- c('c005', 'c010', 'c015', 'c020', 'c025', 'c030', 'c035', 'c040')
ror.med <- c('c045', 'c050', 'c055', 'c060')
ror.high <- c('c065', 'c070', 'c075', 'c080', 'c085', 'c090', 'c095')
patient_annotation_ror <- patient_annotation_ror %>% 
                              mutate(ror.red.hl.class = case_when(ror.red.class %in% ror.low ~ "Low",
                                                                  ror.red.class %in% ror.med ~ "Medium",
                                                                  ror.red.class %in% ror.high ~ "High",
                                                                  TRUE ~ ror.red.class))
patient_annotation_ror <- patient_annotation_ror %>% 
                              mutate(ror.all.hl.class = case_when(ror.all.class %in% ror.low ~ "Low",
                                                                  ror.all.class %in% ror.med ~ "Medium",
                                                                  ror.all.class %in% ror.high ~ "High",
                                                                  TRUE ~ ror.all.class))
patient_annotation_ror$ror.red.hl.class <- factor(patient_annotation_ror$ror.red.hl.class, 
                                                  levels = c("High", "Medium", "Low"))
patient_annotation_ror$ror.all.hl.class <- factor(patient_annotation_ror$ror.all.hl.class, 
                                                  levels = c("High", "Medium", "Low"))

# write.csv(patient_annotation_ror, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ror_all_samples.csv", row.names=FALSE)

patient_annotation_ror_clams <- left_join(patient_annotation_ror, patient_annotation_clams, 
                                          by=c("sample.id", "cancer.type", "dataset"))

# compare classification between predictors
# reduced vs all, all classes
red_vs_all <- table(patient_annotation_ror_clams$ror.red.class, patient_annotation_ror_clams$ror.all.class, useNA='ifany')
red_vs_all
red_vs_all %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()
heatmap(red_vs_all %>% prop.table() %>% '*'(100) %>% round(2), Colv = NA, Rowv = NA,
        col=(viridis(10, direction =-1)))

# reduced vs all, low/medium/high
red_vs_all_hl <- table(patient_annotation_ror_clams$ror.red.hl.class, patient_annotation_ror_clams$ror.all.hl.class, useNA='ifany')
red_vs_all_hl
red_vs_all_hl %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()


# compare classification with CLAMS
# reduced, low/medium/high
clams_vs_red <- table(patient_annotation_ror_clams$clams.class, patient_annotation_ror_clams$ror.red.hl.class, useNA='ifany')
clams_vs_red
clams_vs_red %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# all, low/medium/high
clams_vs_all <- table(patient_annotation_ror_clams$clams.class, patient_annotation_ror_clams$ror.all.hl.class, useNA='ifany')
clams_vs_all
clams_vs_all %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()


# see all classes by CLAMS
patient_annotation_ror_clams %>%
  ggplot(aes(x=ror.red.class, fill=clams.class)) +
  geom_histogram(position = 'identity', stat="count") +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1")) +
  theme_minimal() +
  facet_wrap(~clams.class) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  labs(y = "Number of samples", x = "ROR reduced classification")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ROR_reduced_clams.png", width=9.3, height=5, dpi=300)

patient_annotation_ror_clams %>%
  ggplot(aes(x=ror.all.class, fill=clams.class)) +
  geom_histogram(position = 'identity', stat="count") +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1")) +
  theme_minimal() +
  facet_wrap(~clams.class) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  labs(y = "Number of samples", x = "ROR not reduced classification")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ROR_not_reduced_clams.png", width=9.3, height=5, dpi=300)


# further look at the TRU that are ROR medium risk
# reduced
patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.red.class %in% c("c040", "c045", "c050", "c055")) %>%
  group_by(dataset, cancer.type) %>% tally() # 55/64 in LUAD
tru_medium_red <- patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.red.class %in% c("c040", "c045", "c050", "c055")) %>%
  pull(sample.id)

# all
patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.all.class %in% c("c040", "c045", "c050", "c055")) %>%
  group_by(dataset, cancer.type) %>% tally() # 72/81 in LUAD
tru_medium_all <- patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.all.class %in% c("c040", "c045", "c050", "c055")) %>%
  pull(sample.id)

# common samples
intersect(tru_medium_red, tru_medium_all) #60/507 samples

# check the probabilities of the LUAD ones
# load: ror_results.RData, functions.R
luad_red_prob <- get_samples_prob_from_predictor_result(ror_red_result_tcga_LUAD)
luad_tru_med_red_prob <- subset(luad_red_prob, luad_red_prob$sample.id %in% tru_medium_red)
hist(luad_tru_med_red_prob$sample.prob) # above 0.92
luad_all_prob <- get_samples_prob_from_predictor_result(ror_all_result_tcga_LUAD)
luad_tru_med_all_prob <- subset(luad_all_prob, luad_red_prob$sample.id %in% tru_medium_all)
hist(luad_tru_med_all_prob$sample.prob) # 0.67, 0.79, above 0.95
# probabilities are high, not the problem


# further look at the NonTRU that are super low ROR
# red
patient_annotation_ror_clams %>%
  subset(clams.class == "NonTRU" & ror.red.class == "c005") %>%
  group_by(dataset, cancer.type) %>% tally() %>% arrange(desc(n))

# all
patient_annotation_ror_clams %>%
  subset(clams.class == "NonTRU" & ror.all.class == "c005") %>%
  group_by(dataset, cancer.type) %>% tally() %>% arrange(desc(n))


# BRCA only ------
brca_subset_ror_clams <- subset(patient_annotation_ror_clams, cancer.type == "BRCA")
brca_subset_ror_clams %>% group_by(dataset) %>% tally() #numbers match

# reduced vs all, low/medium/high
red_vs_all_hl_brca <- table(brca_subset_ror_clams$ror.red.hl.class, brca_subset_ror_clams$ror.all.hl.class, useNA='ifany')
red_vs_all_hl_brca
red_vs_all_hl_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# compare classification with CLAMS
# reduced, low/medium/high
clams_vs_red_brca <- table(brca_subset_ror_clams$clams.class, brca_subset_ror_clams$ror.red.hl.class, useNA='ifany')
clams_vs_red_brca
clams_vs_red_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# only TRU
tru_brca <- subset(brca_subset_ror_clams, clams.class == "TRU")
tru_vs_red_brca <- table(tru_brca$clams.class, tru_brca$ror.red.hl.class, useNA='ifany')
tru_vs_red_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# all, low/medium/high
red_vs_all_hl_tru_brca <- table(tru_brca$ror.red.hl.class, tru_brca$ror.all.hl.class, useNA='ifany')
red_vs_all_hl_tru_brca
red_vs_all_hl_tru_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()


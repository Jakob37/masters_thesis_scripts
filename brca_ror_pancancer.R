#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO PREPARE PANCANCER DATA, RUN ROR SSP, ANALYSE DATA
#
#


# CLEAR ENVIRONMENT ----------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(viridisLite)
library(org.Hs.eg.db)

# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION - GOBO
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

# see which genes from predictor are in gex_matrix and filter matrix to keep only those
red_genes <- get.all.pairs.genes(ror.red.aims.gs$all.pairs)
all_genes <- get.all.pairs.genes(ror.all.aims.gs$all.pairs)
red_genes_in_gex <- intersect(red_genes, entrez_ids_gobo)
all_genes_in_gex <- intersect(all_genes, entrez_ids_gobo)
gex_matrix_gobo_red <- gex_matrix_gobo[red_genes_in_gex,]
gex_matrix_gobo_all <- gex_matrix_gobo[all_genes_in_gex,]

# run ROR predictor
ror_red_result_gobo <- applyAIMS(gex_matrix_gobo_red, row.names(gex_matrix_gobo_red), ror.red.aims.gs)
ror_all_result_gobo <- applyAIMS(gex_matrix_gobo_all, row.names(gex_matrix_gobo_all), ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_gobo <- get_samples_class_from_predictor_result(ror_red_result_gobo)
ror_red_class_gobo <- ror_red_class_gobo %>% dplyr::rename(ror.red.class = sample.class)
ror_all_class_gobo <- get_samples_class_from_predictor_result(ror_all_result_gobo)
ror_all_class_gobo <- ror_all_class_gobo %>% dplyr::rename(ror.all.class = sample.class)

# save to patient id table
ror_class_gobo <- inner_join(patient_ids_all_datasets, ror_red_class_gobo, by="sample.id")
ror_class_gobo <- left_join(ror_class_gobo, ror_all_class_gobo, by="sample.id")

rm(gex_matrix_gobo)
rm(entrez_ids_gobo)


# SCAN-B

# Create an entrez id list from the gene symbols present in the gex_matrix  
gene_symbols_scanb <- gex_matrix_scanb[,1:2]
gene_symbols_scanb <- rownames_to_column(data.frame(gene_symbols_scanb), var = "geneSymbol")
gene_symbols_scanb <- gene_symbols_scanb$geneSymbol

entrez_ids_scanb <- mapIds(org.Hs.eg.db, gene_symbols_scanb, 'ENTREZID', 'SYMBOL')
entrez_ids_scanb <- unname(entrez_ids_scanb)
entrez_ids_scanb <- paste0("e",entrez_ids_scanb)

# Rename gex_matrix rows with entrez ids instead of gene symbols
row.names(gex_matrix_scanb) <- entrez_ids_scanb

# see which genes from predictor are in gex_matrix and filter matrix to keep only those
red_genes <- get.all.pairs.genes(ror.red.aims.gs$all.pairs)
all_genes <- get.all.pairs.genes(ror.all.aims.gs$all.pairs)
red_genes_in_gex <- intersect(red_genes, entrez_ids_scanb)
all_genes_in_gex <- intersect(all_genes, entrez_ids_scanb)
gex_matrix_scanb_red <- gex_matrix_scanb[red_genes_in_gex,]
gex_matrix_scanb_all <- gex_matrix_scanb[all_genes_in_gex,]

# run ROR predictor
ror_red_result_scanb <- applyAIMS(gex_matrix_scanb_red, row.names(gex_matrix_scanb_red), ror.red.aims.gs)
ror_all_result_scanb <- applyAIMS(gex_matrix_scanb_all, row.names(gex_matrix_scanb_all), ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_scanb <- get_samples_class_from_predictor_result(ror_red_result_scanb)
ror_red_class_scanb <- ror_red_class_scanb %>% dplyr::rename(ror.red.class = sample.class)
ror_all_class_scanb <- get_samples_class_from_predictor_result(ror_all_result_scanb)
ror_all_class_scanb <- ror_all_class_scanb %>% dplyr::rename(ror.all.class = sample.class)

# save to patient id table
ror_class_scanb <- inner_join(patient_ids_all_datasets, ror_red_class_scanb, by="sample.id")
ror_class_scanb <- left_join(ror_class_scanb, ror_all_class_scanb, by="sample.id")

rm(gex_matrix_scanb)
rm(entrez_ids_scanb)


# TCGA

# Create an entrez id list from the ensembl ids present in the gex_matrix
ensembl_ids_tcga <- gex_matrix_tcga[,1:2]
ensembl_ids_tcga <- rownames_to_column(data.frame(ensembl_ids_tcga), var = "ensemblId")
ensembl_ids_tcga <- ensembl_ids_tcga$ensemblId

entrez_ids_tcga <- mapIds(org.Hs.eg.db, ensembl_ids_tcga, 'ENTREZID', 'ENSEMBL')
entrez_ids_tcga <- unname(entrez_ids_tcga)
entrez_ids_tcga <- paste0("e",entrez_ids_tcga)

# Rename gex_matrix rows with entrez ids instead of reporter ids
row.names(gex_matrix_tcga) <- entrez_ids_tcga

# see which genes from predictor are in gex_matrix and filter matrix to keep only those
red_genes <- get.all.pairs.genes(ror.red.aims.gs$all.pairs)
all_genes <- get.all.pairs.genes(ror.all.aims.gs$all.pairs)
red_genes_in_gex <- intersect(red_genes, entrez_ids_tcga)
all_genes_in_gex <- intersect(all_genes, entrez_ids_tcga)
gex_matrix_tcga_red <- gex_matrix_tcga[red_genes_in_gex,]
gex_matrix_tcga_all <- gex_matrix_tcga[all_genes_in_gex,]

# run ROR predictor for everything
ror_red_result_tcga <- applyAIMS(gex_matrix_tcga_red, row.names(gex_matrix_tcga_red), ror.red.aims.gs)
ror_all_result_tcga <- applyAIMS(gex_matrix_tcga_all, row.names(gex_matrix_tcga_all), ror.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
ror_red_class_tcga <- get_samples_class_from_predictor_result(ror_red_result_tcga)
ror_red_class_tcga <- ror_red_class_tcga %>% dplyr::rename(ror.red.class = sample.class)
ror_all_class_tcga <- get_samples_class_from_predictor_result(ror_all_result_tcga)
ror_all_class_tcga <- ror_all_class_tcga %>% dplyr::rename(ror.all.class = sample.class)

# save to patient id table
ror_class_tcga <- inner_join(patient_ids_all_datasets, ror_red_class_tcga, by="sample.id")
ror_class_tcga <- inner_join(ror_class_tcga, ror_all_class_tcga, by="sample.id")


# run ROR predictor for only BRCA
brca_tcga_subset_samples <- subset(patient_ids_all_datasets, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
col_num <- which(colnames(gex_matrix_tcga) %in% brca_tcga_subset_samples)
brca_tcga_gex_matrix <- gex_matrix_tcga[,c(col_num)]
entrez_ids_tcga_brca <- row.names(brca_tcga_gex_matrix)

ror_red_result_brca_tcga <- applyAIMS(brca_tcga_gex_matrix, entrez_ids_tcga_brca, ror.red.aims.gs)
ror_all_result_brca_tcga <- applyAIMS(brca_tcga_gex_matrix, entrez_ids_tcga_brca, ror.all.aims.gs)

rm(gex_matrix_tcga)
rm(brca_tcga_gex_matrix)
rm(entrez_ids_tcga)


# combine all datasets
patient_annotation_ror <- bind_rows(ror_class_gobo, ror_class_scanb, ror_class_tcga)
patient_annotation_ror$ror.red.class <- factor(patient_annotation_ror$ror.red.class)
patient_annotation_ror$ror.all.class <- factor(patient_annotation_ror$ror.all.class)


# compare what is low (<40) and what is High (>60)
ror.low <- c('c005', 'c010', 'c015', 'c020', 'c025', 'c030', 'c035', 'c040')
ror.med <- c('c045', 'c050', 'c055', 'c060')
ror.high <- c('c065', 'c070', 'c075', 'c080', 'c085', 'c090', 'c095')
patient_annotation_ror$ror.red.class <- as.character(patient_annotation_ror$ror.red.class)
patient_annotation_ror <- patient_annotation_ror %>% 
  mutate(ror.red.hl.class = case_when(ror.red.class %in% ror.low ~ "Low",
                                      ror.red.class %in% ror.med ~ "Medium",
                                      ror.red.class %in% ror.high ~ "High",
                                      TRUE ~ ror.red.class))
patient_annotation_ror$ror.all.class <- as.character(patient_annotation_ror$ror.all.class)
patient_annotation_ror <- patient_annotation_ror %>% 
  mutate(ror.all.hl.class = case_when(ror.all.class %in% ror.low ~ "Low",
                                      ror.all.class %in% ror.med ~ "Medium",
                                      ror.all.class %in% ror.high ~ "High",
                                      TRUE ~ ror.all.class))

# RESULT
# Save final table
# write.csv(patient_annotation_ror, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/ror_all_samples.csv", row.names=FALSE)


setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/")
#save(ror_red_result_gobo, ror_all_result_gobo, ror_red_result_scanb, ror_all_result_scanb, 
#      ror_red_result_tcga, ror_all_result_tcga, ror_red_result_brca_tcga, ror_all_result_brca_tcga, file='ror_results.RData')


# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R
# ROR RESULTS 
# CLAMS RESULTS


# Pancancer  -------

# ROR by cancer type
patient_annotation_ror <- patient_annotation_ror %>% mutate(group.to.analyze = paste(cancer.type, dataset, sep="_"))
table(patient_annotation_ror$group.to.analyze, patient_annotation_ror$ror.red.hl.class) %>% addmargins() %>% print(zero.print=".")
#              High   Low Medium   Sum
# ACC_TCGA       11    34     33    78
# BLCA_TCGA     336    10     58   404
# BRCA_GOBO     264   524   1093  1881
# BRCA_SCAN-B  1004  1669    847  3520
# BRCA_TCGA     655   161    256  1072
# CESC_TCGA     299     .      2   301
# CHOL_TCGA       7     7     22    36
# COAD_TCGA     252     .     24   276
# DLBC_TCGA      46     .      2    48
# ESCA_TCGA     135     2     22   159
# GBM_TCGA       26    39     87   152
# HNSC_TCGA     456     2     36   494
# KICH_TCGA       2    57      5    64
# KIRC_TCGA      13   439     60   512
# KIRP_TCGA       3   243     37   283
# LGG_TCGA       13   457     38   508
# LIHC_TCGA      93    92    183   368
# LUAD_TCGA     242   100    165   507
# LUSC_TCGA     432     5     49   486
# MESO_TCGA      30    11     45    86
# OV_TCGA       175     3     56   234
# PAAD_TCGA      21    49     85   155
# PCPG_TCGA       .   177      .   177
# PRAD_TCGA       6   433     52   491
# READ_TCGA      77     .     12    89
# SARC_TCGA     115    48     91   254
# SKCM_TCGA      63     2     38   103
# STAD_TCGA     260    30     72   362
# TGCT_TCGA     128     2      3   133
# THCA_TCGA       .   489      8   497
# THYM_TCGA      76    20     23   119
# UCEC_TCGA     146     7     18   171
# UCS_TCGA       52     .      4    56
# UVM_TCGA        .    36     44    80
# Sum          5438  5148   3570 14156

patient_annotation_ror_clams <- left_join(patient_annotation_ror, patient_annotation_clams, 
                                          by=c("sample.id", "cancer.type", "dataset"))
patient_annotation_ror_clams$ror.red.hl.class <- factor(patient_annotation_ror_clams$ror.red.hl.class, 
                                                  levels = c("High", "Medium", "Low"))
patient_annotation_ror_clams$ror.all.hl.class <- factor(patient_annotation_ror_clams$ror.all.hl.class, 
                                                  levels = c("High", "Medium", "Low"))


# compare classification between predictors, pancancer
# reduced vs all, all c0 classes
red_vs_all <- table(patient_annotation_ror_clams$ror.red.class, patient_annotation_ror_clams$ror.all.class, useNA='ifany')
red_vs_all %>% print(zero.print=".")
red_vs_all %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins() %>% print(zero.print=".")
heatmap(red_vs_all, Colv = NA, Rowv = NA, col=(viridis(10, direction =-1)))

# reduced vs all, low/medium/high
red_vs_all_hl <- table(patient_annotation_ror_clams$ror.red.hl.class, patient_annotation_ror_clams$ror.all.hl.class, useNA='ifany')
red_vs_all_hl %>% print(zero.print=".")
#        High Medium  Low
# High   5381     57    .
# Medium  812   2579  179
# Low       .    174 4974
# 12934/14156 (91.37%) are the same between both methods
red_vs_all_hl %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# compare classification with CLAMS, pancancer
# reduced, low/medium/high
clams_vs_red <- table(patient_annotation_ror_clams$clams.class, patient_annotation_ror_clams$ror.red.hl.class, useNA='ifany')
clams_vs_red %>% print(zero.print=".")
#         High Medium  Low
# NonTRU 5438   3509 3358
# TRU       .     61 1790
# 1790/1851 (96.7%) of TRU samples are Low
clams_vs_red %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# all, low/medium/high
clams_vs_all <- table(patient_annotation_ror_clams$clams.class, patient_annotation_ror_clams$ror.all.hl.class, useNA='ifany')
clams_vs_all %>% print(zero.print=".")
#        High Medium  Low
# NonTRU 6192   2762 3351
# TRU       1     48 1802
# 1802/1851 (97.35%) of TRU samples are Low
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


# further look at the TRU that are ROR "medium" risk (c040 is low)
# reduced
patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.red.class %in% c("c040", "c045", "c050", "c055", "c060")) %>%
  group_by(dataset, cancer.type) %>% tally() # 56/65 in LUAD
# dataset cancer.type       n
# 1 GOBO    BRCA            4
# 2 TCGA    LUAD           56
# 3 TCGA    LUSC            3
# 4 TCGA    THCA            2
tru_40to60_red <- patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.red.class %in% c("c040", "c045", "c050", "c055", "c060")) %>%
  pull(sample.id)

# all
patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.all.class %in% c("c040", "c045", "c050", "c055", "c060")) %>%
  group_by(dataset, cancer.type) %>% tally() # 71/80 in LUAD
# dataset cancer.type     n
# 1 GOBO    BRCA            2
# 2 TCGA    KIRP            1
# 3 TCGA    LIHC            2
# 4 TCGA    LUAD           71
# 5 TCGA    LUSC            2
# 6 TCGA    THCA            2

tru_40to60_all <- patient_annotation_ror_clams %>%
  subset(clams.class == "TRU" & ror.all.class %in% c("c040", "c045", "c050", "c055", "c060")) %>%
  pull(sample.id)

# common samples
intersect(tru_40to60_red, tru_40to60_all) #61 samples

# check the probabilities of the LUAD ones - didn't redo it
# load: ror_results.RData, functions.R
luad_red_prob <- get_samples_prob_from_predictor_result(ror_red_result_tcga_LUAD)
luad_tru_4060_red_prob <- subset(luad_red_prob, luad_red_prob$sample.id %in% tru_40to60_red)
hist(luad_tru_4060_red_prob$sample.prob) # above 0.92
luad_all_prob <- get_samples_prob_from_predictor_result(ror_all_result_tcga_LUAD)
luad_tru_4060_all_prob <- subset(luad_all_prob, luad_red_prob$sample.id %in% tru_40to60_all)
hist(luad_tru_4060_all_prob$sample.prob) # 0.67, 0.79, above 0.95
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
red_vs_all_hl_brca %>% print(zero.print=".")
#        High Medium  Low
# High   1893     30    .
# Medium  478   1582  136
# Low       .     74 2280
# so 5755/6473 (88.91%) samples are the same category with both methods
red_vs_all_hl_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# per dataset
# GOBO
brca_gobo_ror <- brca_subset_ror_clams %>% filter(dataset == "GOBO")
table(brca_gobo_ror$ror.red.hl.class, brca_gobo_ror$ror.all.hl.class, useNA='ifany')  %>% print(zero.print=".")
#         High Medium Low
# High    262      2   .
# Medium  270    764  59
# Low       .     25 499
# so 1525/1881 (81.07%) samples are the same category with both methods

# SCAN-B
brca_scanb_ror <- brca_subset_ror_clams %>% filter(dataset == "SCAN-B")
table(brca_scanb_ror$ror.red.hl.class, brca_scanb_ror$ror.all.hl.class, useNA='ifany') %>% print(zero.print=".")
#        High Medium  Low
# High    996      8    .
# Medium  178    600   69
# Low       .     38 1631
# so 3227/3520 (91.68%) samples are the same category with both methods

# TCGA
brca_tcga_ror <- brca_subset_ror_clams %>% filter(dataset == "TCGA")
table(brca_tcga_ror$ror.red.hl.class, brca_tcga_ror$ror.all.hl.class, useNA='ifany') %>% print(zero.print=".")
# red as rows, all as columns
#         High Medium Low
# High    635     20   .
# Medium   30    218   8
# Low       .     11 150
# so 1003/1072 (93.56%) samples are the same category with both methods



# compare classification with CLAMS
# reduced, low/medium/high
clams_vs_red_brca <- table(brca_subset_ror_clams$clams.class, brca_subset_ror_clams$ror.red.hl.class, useNA='ifany')
clams_vs_red_brca %>% print(zero.print=".")
#        High Medium  Low
# NonTRU 1923   2194 1593
# TRU       .      2  761 # 99.74% of BRCA TRU samples are low ROR
clams_vs_red_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()

# per dataset
# GOBO
clams_vs_red_brca_gobo <- table(brca_gobo_ror$clams.class, brca_gobo_ror$ror.red.hl.class, useNA='ifany')
clams_vs_red_brca_gobo %>% print(zero.print=".")
#         High Medium  Low
# NonTRU  264   1091  326
# TRU       .      2  198 # 99% of TRU are low ROR
# SCAN-B
clams_vs_red_brca_scanb <- table(brca_scanb_ror$clams.class, brca_scanb_ror$ror.red.hl.class, useNA='ifany')
clams_vs_red_brca_scanb %>% print(zero.print=".")
#        High Medium  Low
# NonTRU 1004    847 1145
# TRU       .      .  524 # 100% of TRU are low ROR
# TCGA
clams_vs_red_brca_tcga <- table(brca_tcga_ror$clams.class, brca_tcga_ror$ror.red.hl.class, useNA='ifany')
clams_vs_red_brca_tcga %>% print(zero.print=".")
#         High Medium Low
# NonTRU  655    256 122
# TRU       .      .  39 # 100% of TRU are low ROR

# with PROLIFERATION, pancancer

patient_annotation_ror_clams_prolif <- left_join(patient_annotation_ror_clams,
                                                 patient_annotation_prolif,
                                                 by = c("sample.id", "cancer.type", "dataset"))
ror.names <- c("High ROR", "Medium ROR", "Low ROR")
names(ror.names) <- c("High", "Medium", "Low")

patient_annotation_ror_clams_prolif %>%
  # filter(ror.red.hl.class == "Low") %>%
  ggplot(aes(x=clams.class, y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  facet_wrap(~ror.red.hl.class, labeller = labeller(ror.red.hl.class = ror.names)) +
  ggtitle("All datasets")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams.png", width=5, height=4, dpi=300)


# with PROLIFERATION, BRCA only, exclude GOBO

brca_ror_prolif_nogobo <- patient_annotation_ror_clams_prolif %>% filter(cancer.type == "BRCA" & dataset != "GOBO")

sample_count <- group_by(brca_ror_prolif_nogobo, ror.red.hl.class, clams.class) %>% tally()
# ror.red.hl.class clams.class     n
# High             NonTRU       1659
# Low              NonTRU       1267
# Low              TRU           563
# Medium           NonTRU       1103

brca_ror_prolif_nogobo %>%
  # filter(ror.red.hl.class == "Low") %>%
  ggplot(aes(x=clams.class, y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  facet_wrap(~ror.red.hl.class, labeller = labeller(ror.red.hl.class = ror.names)) +
  ggtitle("BRCA only, no GOBO, reduced SSP")
    
# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams_BRCA_nogobo.png", width=5, height=4, dpi=300)
# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams_BRCA_nogobo.pdf", width=5, height=4)

  
  

# per dataset
# GOBO
brca_gobo_ror_prolif <- patient_annotation_ror_clams_prolif %>% filter(dataset == "GOBO")
brca_gobo_ror_prolif %>%
  # filter(ror.red.hl.class == "Low") %>%
  ggplot(aes(x=clams.class, y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  facet_wrap(~ror.red.hl.class, labeller = labeller(ror.red.hl.class = ror.names)) +
  ggtitle("BRCA (GOBO)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams_BRCA_GOBO.png", width=5, height=4, dpi=300)

# SCAN-B
brca_scanb_ror_prolif <- patient_annotation_ror_clams_prolif %>% filter(dataset == "SCAN-B")
brca_scanb_ror_prolif %>%
  # filter(ror.red.hl.class == "Low") %>%
  ggplot(aes(x=clams.class, y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  facet_wrap(~ror.red.hl.class, labeller = labeller(ror.red.hl.class = ror.names)) +
  ggtitle("BRCA (SCAN-B)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams_BRCA_SCANB.png", width=5, height=4, dpi=300)

# TCGA
brca_tcga_ror_prolif <- patient_annotation_ror_clams_prolif %>% filter(dataset == "TCGA" & cancer.type == "BRCA")
brca_tcga_ror_prolif %>%
  # filter(ror.red.hl.class == "Low") %>%
  ggplot(aes(x=clams.class, y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  facet_wrap(~ror.red.hl.class, labeller = labeller(ror.red.hl.class = ror.names)) +
  ggtitle("BRCA (TCGA)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/prolif/ROR_red_clams_BRCA_TCGA.png", width=5, height=4, dpi=300)
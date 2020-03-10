#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO PREPARE BREAST CANCER DATA, RUN PAM50 SSP, ANALYSE DATA
#
#


# CLEAR ENVIRONMENT ----------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(org.Hs.eg.db)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX
# FUNCTION TO RUN DIFFERENT SSPS run_ssps_function.R
# functions.R
# PAM50 SSPs

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/")


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

# run PAM50 predictor
pam50_red_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, pam50.red.aims.gs)
pam50_all_result_gobo <- applyAIMS(gex_matrix_gobo, entrez_ids_gobo, pam50.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
pam50_red_class_gobo <- get_samples_class_from_predictor_result(pam50_red_result_gobo)
pam50_red_class_gobo <- pam50_red_class_gobo %>% dplyr::rename(pam50.red.class = sample.class)
pam50_all_class_gobo <- get_samples_class_from_predictor_result(pam50_all_result_gobo)
pam50_all_class_gobo <- pam50_all_class_gobo %>% dplyr::rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_gobo <- inner_join(patient_ids_all_datasets, pam50_red_class_gobo, by="sample.id")
pam50_class_gobo <- inner_join(pam50_class_gobo, pam50_all_class_gobo, by="sample.id")

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

# run PAM50 predictor
pam50_red_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, pam50.red.aims.gs)
pam50_all_result_scanb <- applyAIMS(gex_matrix_scanb, entrez_ids_scanb, pam50.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
pam50_red_class_scanb <- get_samples_class_from_predictor_result(pam50_red_result_scanb)
pam50_red_class_scanb <- pam50_red_class_scanb %>% dplyr::rename(pam50.red.class = sample.class)
pam50_all_class_scanb <- get_samples_class_from_predictor_result(pam50_all_result_scanb)
pam50_all_class_scanb <- pam50_all_class_scanb %>% dplyr::rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_scanb <- inner_join(patient_ids_all_datasets, pam50_red_class_scanb, by="sample.id")
pam50_class_scanb <- inner_join(pam50_class_scanb, pam50_all_class_scanb, by="sample.id")

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

# Filter matrix to keep only genes the SSP will use
red_genes <- get.unique.genes(pam50.red.aims.gs[["all.pairs"]])
all_genes <- get.unique.genes(pam50.all.aims.gs[["all.pairs"]])
red_all_genes <- union(red_genes, all_genes)
gex_matrix_tcga <- gex_matrix_tcga[red_all_genes,]

# run PAM50 predictor for everything
pam50_red_result_tcga <- applyAIMS(gex_matrix_tcga, red_all_genes, pam50.red.aims.gs)
pam50_all_result_tcga <- applyAIMS(gex_matrix_tcga, red_all_genes, pam50.all.aims.gs)

# extract the classification per sample from the predictor result, rename column to reflect aims.gs
pam50_red_class_tcga <- get_samples_class_from_predictor_result(pam50_red_result_tcga)
pam50_red_class_tcga <- pam50_red_class_tcga %>% dplyr::rename(pam50.red.class = sample.class)
pam50_all_class_tcga <- get_samples_class_from_predictor_result(pam50_all_result_tcga)
pam50_all_class_tcga <- pam50_all_class_tcga %>% dplyr::rename(pam50.all.class = sample.class)

# save to patient id table
pam50_class_tcga <- inner_join(patient_ids_all_datasets, pam50_red_class_tcga, by="sample.id")
pam50_class_tcga <- inner_join(pam50_class_tcga, pam50_all_class_tcga, by="sample.id")

# run ROR predictor for only BRCA
brca_tcga_subset_samples <- subset(patient_ids_all_datasets, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
col_num <- which(colnames(gex_matrix_tcga) %in% brca_tcga_subset_samples)
brca_tcga_gex_matrix <- gex_matrix_tcga[,c(col_num)]
entrez_ids_tcga_brca <- row.names(brca_tcga_gex_matrix)

pam50_red_result_brca_tcga <- applyAIMS(brca_tcga_gex_matrix, entrez_ids_tcga_brca, pam50.red.aims.gs)
pam50_all_result_brca_tcga <- applyAIMS(brca_tcga_gex_matrix, entrez_ids_tcga_brca, pam50.all.aims.gs)

rm(gex_matrix_tcga)
rm(entrez_ids_tcga)

# combine all datasets
patient_annotation_pam50 <- bind_rows(pam50_class_gobo, pam50_class_scanb, pam50_class_tcga)
patient_annotation_pam50$pam50.red.class <- factor(patient_annotation_pam50$pam50.red.class)
patient_annotation_pam50$pam50.all.class <- factor(patient_annotation_pam50$pam50.all.class)

# RESULT
# Save final table
write.csv(patient_annotation_pam50, "/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/pam50_all_samples.csv", row.names=FALSE)

# setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/")
# save(pam50_red_result_gobo, pam50_all_result_gobo, pam50_red_result_scanb, pam50_all_result_scanb,
#      pam50_red_result_tcga, pam50_all_result_tcga, pam50_red_result_brca_tcga, pam50_all_result_brca_tcga, file='pam50_results.RData')


# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# PAM50 RESULTS 
# CLAMS RESULTS

# Pancancer  -------
# compare classification between predictors
table(patient_annotation_pam50$pam50.red.class, 
        patient_annotation_pam50$pam50.all.class, useNA='ifany') %>% print(zero.print=".")
#        Basal Her2 LumA LumB
# Basal  3158   64    .    .
# Her2   1756 4430   84  310
# LumA      4  115 2539  186
# LumB      .  103   37 1370
# 11497/14156 (81.22%) samples are classified the same between both predictors
table(patient_annotation_pam50$pam50.red.class, 
      patient_annotation_pam50$pam50.all.class, useNA='ifany') %>% 
      prop.table() %>% '*'(100) %>% round(2) %>% addmargins() %>% print(zero.print=".")
#       Basal  Her2  LumA  LumB   Sum
# Basal 22.31  0.45     .     . 22.76
# Her2  12.40 31.29  0.59  2.19 46.47
# LumA   0.03  0.81 17.94  1.31 20.09
# LumB      .  0.73  0.26  9.68 10.67
# Sum   34.74 33.28 18.79 13.18 99.99

# with CLAMS classification as well
# reduced
patient_annotation_pam50_clams <- left_join(patient_annotation_pam50, patient_annotation_clams, 
                                          by=c("sample.id", "cancer.type", "dataset"))
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.red.class, useNA='ifany')
#         Basal Her2 LumA LumB
# NonTRU  3136 5588 2072 1509
# TRU       86  992  772    1
# of 1851 TRU samples, 4.65% are Basal, 53.59% are Her2, 41.71% are LumA, 0.05% are LumB
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.red.class, 
      useNA='ifany') %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()
#         Basal   Her2   LumA   LumB    Sum
# NonTRU  22.15  39.47  14.64  10.66  86.92
# TRU      0.61   7.01   5.45   0.01  13.08
# Sum     22.76  46.48  20.09  10.67 100.00

# all 19k
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.all.class, useNA='ifany')
#         Basal Her2 LumA LumB
# NonTRU  4468 4083 1891 1863
# TRU      450  629  769    3
# of 1851 TRU samples, 24.31% are Basal, 33.98% are Her2, 41.55% are LumA, 0.16% are LumB
table(patient_annotation_pam50_clams$clams.class, patient_annotation_pam50_clams$pam50.all.class, 
      useNA='ifany') %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()
#        Basal  Her2  LumA  LumB   Sum
# NonTRU 31.56 28.84 13.36 13.16 86.92
# TRU     3.18  4.44  5.43  0.02 13.07
# Sum    34.74 33.28 18.79 13.18 99.99

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
  scale_y_continuous(breaks = c(1000,2000,3000,4000,5000,6000), limits = c(0,6000)) +
  labs(y = "Number of samples", x = "PAM50 reduced classification")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/PAM50_reduced_clams.png", width=9.3, height=5, dpi=300)

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
  scale_y_continuous(breaks = c(1000,2000,3000,4000,5000,6000), limits = c(0,6000)) +
  labs(y = "Number of samples", x = "PAM50 not reduced classification")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/PAM50_not_reduced_clams.png", width=9.3, height=5, dpi=300)


# BRCA only ------
brca_subset_pam50_clams <- subset(patient_annotation_pam50_clams, cancer.type == "BRCA")
brca_subset_pam50_clams %>% group_by(dataset) %>% tally() #numbers match

# reduced vs all, all datasets together
red_vs_all_brca <- table(brca_subset_pam50_clams$pam50.red.class, 
                         brca_subset_pam50_clams$pam50.all.class, useNA='ifany') %>% print(zero.print=".")
print(red_vs_all_brca %>% addmargins(), zero.print = ".")
#       Basal Her2 LumA LumB  Sum
# Basal   761    4    .    .  765
# Her2     68 1051   58  310 1487
# LumA      4   39 2483  186 2712
# LumB      .  102   37 1370 1509
# Sum     833 1196 2578 1866 6473
# 5665/6473 (87.52%) BRCA samples are classified as the same between predictors
print(red_vs_all_brca %>% prop.table() %>% '*'(100) %>% round(1) %>% addmargins(), zero.print = ".")
#       Basal  Her2  LumA  LumB   Sum
# Basal  11.8   0.1     .     .  11.9
# Her2    1.1  16.2   0.9   4.8  23.0
# LumA    0.1   0.6  38.4   2.9  42.0
# LumB      .   1.6   0.6  21.2  23.4
# Sum    13.0  18.5  39.9  28.9 100.3

# reduced vs all, per dataset
datasets <- levels(factor(brca_subset_pam50_clams$dataset))
for (my_dataset in datasets) {
  print(my_dataset)
  subset_dataset <- subset(brca_subset_pam50_clams, dataset == my_dataset)
  print(table(subset_dataset$pam50.red.class, subset_dataset$pam50.all.class, useNA='ifany') %>% addmargins(), zero.print = ".")
}
# GOBO
#       Basal Her2 LumA LumB  Sum
# Basal   223    .    .    .  223
# Her2     48  422   46  275  791
# LumA      .   13  454  148  615
# LumB      .   10    5  237  252
# Sum     271  445  505  660 1881
# 1336/1881 (71.03%) are concordant

# SCAN-B
#       Basal Her2 LumA LumB  Sum
# Basal   351    4    .    .  355
# Her2     15  434    8   15  472
# LumA      4   23 1784   31 1842
# LumB      .   53   30  768  851
# Sum     370  514 1822  814 3520
# 3337/3520 (94.8%) are concordant

# TCGA
#       Basal Her2 LumA LumB  Sum
# Basal   187    .    .    .  187
# Her2      5  195    4   20  224
# LumA      .    3  245    7  255
# LumB      .   39    2  365  406
# Sum     192  237  251  392 1072
# 992/1072 (92.54%) are concordant

# compare classification with CLAMS
# reduced
clams_vs_red_brca <- table(brca_subset_pam50_clams$clams.class, brca_subset_pam50_clams$pam50.red.class, useNA='ifany')
clams_vs_red_brca
#         Basal Her2 LumA LumB
# NonTRU   757 1449 1996 1508
# TRU        8   38  716    1
# of 763 BRCA TRU samples, 1.05% are Basal, 4.98% are Her2, 93.84% are LumA, 0.13% are LumB
clams_vs_red_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()
#         Basal   Her2   LumA   LumB    Sum
# NonTRU  11.69  22.39  30.84  23.30  88.22
# TRU      0.12   0.59  11.06   0.02  11.79
# Sum     11.81  22.98  41.90  23.32 100.01

# all
clams_vs_all_brca <- table(brca_subset_pam50_clams$clams.class, brca_subset_pam50_clams$pam50.all.class, useNA='ifany')
clams_vs_all_brca
#         Basal Her2 LumA LumB
# NonTRU   813 1179 1855 1863
# TRU       20   17  723    3
# of 763 BRCA TRU samples, 2.62% are Basal, 2.23% are Her2, 94.76% are LumA, 0.39% are LumB
clams_vs_all_brca %>% prop.table() %>% '*'(100) %>% round(2) %>% addmargins()
#         Basal   Her2   LumA   LumB    Sum
# NonTRU  12.56  18.21  28.66  28.78  88.21
# TRU      0.31   0.26  11.17   0.05  11.79
# Sum     12.87  18.47  39.83  28.83 100.00


# with PROLIFERATION

patient_annotation_pam50_clams_prolif <- left_join(patient_annotation_pam50_clams,
                                                 patient_annotation_prolif,
                                                 by = c("sample.id", "cancer.type", "dataset"))
patient_annotation_pam50_clams_prolif$pam50.red.class <- factor(patient_annotation_pam50_clams_prolif$pam50.red.class,
                                                                 levels=c("Basal", "Her2", "LumB", "LumA"))
pam.names <- c("Basal", "HER2", "LumB", "LumA")
names(pam.names) <- c("Basal", "Her2", "LumB", "LumA")

# pancancer
patient_annotation_pam50_clams_prolif %>%
  # filter(pam50.red.class == "LumA") %>%
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
  facet_wrap(~pam50.red.class, labeller = labeller(pam50.red.class = pam.names), ncol=4) +
  ggtitle("Pancancer")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/prolif/PAM50_red_clams.png", width=6, height=3.5, dpi=300)


# BRCA only
brca_pam50_prolif <- patient_annotation_pam50_clams_prolif %>% filter(cancer.type == "BRCA")

brca_pam50_prolif %>%
  # filter(pam50.red.class == "LumA") %>%
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
  facet_wrap(~pam50.red.class, labeller = labeller(pam50.red.class = pam.names), ncol=4) +
  ggtitle("BRCA only, all datasets")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/prolif/PAM50_brca_red_clams.png", width=6, height=3.5, dpi=300)

# by dataset
brca_gobo_pam50_prolif <- brca_pam50_prolif %>% filter(dataset == "GOBO")
brca_gobo_pam50_prolif %>%
  # filter(pam50.red.class == "LumA") %>%
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
  facet_wrap(~pam50.red.class, labeller = labeller(pam50.red.class = pam.names), ncol=4) +
  ggtitle("BRCA (GOBO)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/prolif/PAM50_brca_gobo_red_clams.png", width=6, height=3.5, dpi=300)

brca_scanb_pam50_prolif <- brca_pam50_prolif %>% filter(dataset == "SCAN-B")
brca_scanb_pam50_prolif %>%
  # filter(pam50.red.class == "LumA") %>%
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
  facet_wrap(~pam50.red.class, labeller = labeller(pam50.red.class = pam.names), ncol=4) +
  ggtitle("BRCA (SCAN-B)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/prolif/PAM50_brca_scanb_red_clams.png", width=6, height=3.5, dpi=300)

brca_tcga_pam50_prolif <- brca_pam50_prolif %>% filter(dataset == "TCGA" & cancer.type == "BRCA")
brca_tcga_pam50_prolif %>%
  # filter(pam50.red.class == "LumA") %>%
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
  facet_wrap(~pam50.red.class, labeller = labeller(pam50.red.class = pam.names), ncol=4) +
  ggtitle("BRCA (TCGA)")

# ggsave("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/prolif/PAM50_brca_tcga_red_clams.png", width=6, height=3.5, dpi=300)


# Comparing the probabilities  ------

# PAM 50, reduced by dataset
# load: pam50_results.RData, functions.R
# GOBO
red_prob_gobo <- get_samples_prob_from_predictor_result(pam50_red_result_gobo)
hist(red_prob_gobo$sample.prob)
mean(red_prob_gobo$sample.prob) # 0.9933015

# SCAN-B
red_prob_scanb <- get_samples_prob_from_predictor_result(pam50_red_result_scanb)
hist(red_prob_scanb$sample.prob)
mean(red_prob_scanb$sample.prob) # 0.9966764

# TCGA
red_prob_tcga <- get_samples_prob_from_predictor_result(pam50_red_result_tcga)
hist(red_prob_tcga$sample.prob)
mean(red_prob_tcga$sample.prob) # 0.9871545

# PAM 50, NOT reduced by dataset
# GOBO
all_prob_gobo <- get_samples_prob_from_predictor_result(pam50_all_result_gobo)
hist(all_prob_gobo$sample.prob)
mean(all_prob_gobo$sample.prob) # 0.9875419

# SCAN-B
all_prob_scanb <- get_samples_prob_from_predictor_result(pam50_all_result_scanb)
hist(all_prob_scanb$sample.prob)
mean(all_prob_scanb$sample.prob) # 0.995342

# TCGA
all_prob_tcga <- get_samples_prob_from_predictor_result(pam50_all_result_tcga)
hist(all_prob_tcga$sample.prob)
mean(all_prob_tcga$sample.prob) # 0.9908023

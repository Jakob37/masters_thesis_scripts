#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# SCRIPT TO PREPARE PAN CANCER DATA, GET IMMUNE SIGNATURE MEASURE AND PLOT
#
#

# CLEAR ENVIRONMENT ------------------------------------------------------------------------------------------------------------------------
  
rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ----------------------------------------------------------------------------------------------------------------

library(tidyverse)



# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1) PATIENT IDS FROM ANNOTATION
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX


# ANALYSIS  ------------------------------------------------------

setwd("/media/deboraholi/Data/LUND/9 THESIS/2_immune/")


# GET GENE SYMBOLS AND LIST  -------------------------------------

# GOBO
reporter_to_symbol_gobo <- gex_matrix_gobo[,1:2]
reporter_to_symbol_gobo <- rownames_to_column(data.frame(reporter_to_symbol_gobo), var = "reporterId")
reporter_to_symbol_gobo <- left_join(reporter_to_symbol_gobo, gene_table_gobo[,c("reporterId","geneSymbol")])
gene_symbol_column_gobo <- reporter_to_symbol_gobo["geneSymbol"]

rm(gene_table_gobo)

length(gene_symbol_column_gobo$geneSymbol)
length(unique(gene_symbol_column_gobo$geneSymbol))

genes_gobo <- unique(gene_symbol_column_gobo$geneSymbol)

# SCAN-B
# row.names already are the gene symbols, no duplicate genes
reporter_to_symbol_scanb <- gex_matrix_scanb[,1:2]
reporter_to_symbol_scanb <- rownames_to_column(data.frame(gex_matrix_scanb), var = "geneSymbol")
gene_symbol_column_scanb <- reporter_to_symbol_scanb["geneSymbol"]

length(gene_symbol_column_scanb$geneSymbol)
length(unique(gene_symbol_column_scanb$geneSymbol))

genes_scanb <- unique(gene_symbol_column_scanb$geneSymbol)


#TCGA
reporter_to_symbol_tcga <- gex_matrix_tcga[,1:2]
reporter_to_symbol_tcga <- rownames_to_column(data.frame(reporter_to_symbol_tcga), var = "reporterId")
reporter_to_symbol_tcga <- left_join(reporter_to_symbol_tcga, gene_table_tcga[,c("ENSG", "SYMBOL")], by=c("reporterId" = "ENSG"))
reporter_to_symbol_tcga <- reporter_to_symbol_tcga %>% rename(geneSymbol=SYMBOL)
gene_symbol_column_tcga <- reporter_to_symbol_tcga["geneSymbol"]

rm(gene_table_tcga)

length(gene_symbol_column_tcga$geneSymbol)
length(unique(gene_symbol_column_tcga$geneSymbol))

genes_tcga <- unique(gene_symbol_column_tcga$geneSymbol)



# GET COMMON GENE SYMBOLS  ---------------------------------

common_genes <- intersect(intersect(genes_gobo, genes_scanb), genes_tcga) #10,916

length(common_genes)
length(unique(common_genes))

rm(genes_gobo)
rm(genes_scanb)
rm(genes_tcga)



# INPUT MODULE ------------------------------------------

mod_immune <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/2_immune/module/immune_signature.csv")
mod_immune <- mod_immune$geneSymbol
length(mod_immune)
mod_immune

mod_immune_in_common <- mod_immune[mod_immune %in% common_genes]
length(mod_immune_in_common)
mod_immune_in_common

rm(mod_immune)

# Result
  write.table(mod_immune_in_common, file='immune_genes_in_common.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)


# MODULE ANALYSIS  ------------------------------------

# GET SIGNATURE FUNCTION FROM FILE
source("/media/deboraholi/Data/LUND/9 THESIS/src/functions.R")


# GOBO
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_gobo))) {
  sample_column <- gex_matrix_gobo[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_gobo, sample_column, mod_immune_in_common, max)
  patient_ids_all_datasets$immune.value.max[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  rm(sample_column)
  if (i %% 100 == 0) { print(i) }
}

# with sum as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_gobo))) {
  sample_column <- gex_matrix_gobo[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_gobo, sample_column, mod_immune_in_common, sum)
  patient_ids_all_datasets$immune.value.sum[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  if (i %% 100 == 0) { print(i) }
}

rm(gex_matrix_gobo)

# SCAN-B
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_scanb))) {
  sample_column <- gex_matrix_scanb[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_scanb, sample_column, mod_immune_in_common, max)
  patient_ids_all_datasets$immune.value.max[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  if (i %% 100 == 0) { print(i) }
}

# with sum as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_scanb))) {
  sample_column <- gex_matrix_scanb[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_scanb, sample_column, mod_immune_in_common, sum)
  patient_ids_all_datasets$immune.value.sum[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  if (i %% 100 == 0) { print(i) }
}

rm(gex_matrix_scanb)


# TCGA
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_tcga))) {
  sample_column <- gex_matrix_tcga[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_tcga, sample_column, mod_immune_in_common, max)
  patient_ids_all_datasets$immune.value.max[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  if (i %% 200 == 0) { print(i) }
}

# with sum as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_tcga))) {
  sample_column <- gex_matrix_tcga[,i,drop=FALSE]
  immune_value <- get_signature_value(gene_symbol_column_tcga, sample_column, mod_immune_in_common, sum)
  patient_ids_all_datasets$immune.value.sum[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- immune_value
  if (i %% 200 == 0) { print(i) }
}

rm(gex_matrix_tcga)


# RESULT
# Save final table
write.csv(patient_ids_all_datasets, "immune_all_samples.csv", row.names=FALSE)


# LOAD RESULTS FROM data_input.R  -----------------------------------------------

## from here blah --------------------------------------------------------------------------------------------

# CLAMS RESULTS
# IMMUNE ANALYSIS RESULTS



# MERGE RESULTS WITH SPECIFIC COLUMNS FROM PATIENT ANNOTATION  --------------------------------





# PLOTS  ---------------------------------------------------

# GET NAMES FROM FILE
source("/media/deboraholi/Data/LUND/9 THESIS/src/plots_names.R")

# # fred box
# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, fred.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, fred.value), y=fred.value, fill=clams.class)) +
#   geom_boxplot(outlier.size = 0.1) +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   theme(axis.text.x = element_text(angle = 60)) +
#   scale_y_continuous(breaks = c(50000, 150000), labels = c('Lower', 'Higher'))
# 
# # fred box flip
# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, fred.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, fred.value), y=fred.value, fill=clams.class)) +
#   geom_boxplot(outlier.size = 0.1) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   scale_y_continuous(breaks = c(50000, 150000), labels = c('Lower', 'Higher'))
# 
# # fred violin
# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, fred.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, fred.value), y=fred.value, fill=clams.class)) +
#   geom_violin() +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   scale_y_continuous(breaks = c(50000, 150000), labels = c('Lower', 'Higher'))


# karl box - no CLAMS separation
patient_annotation_all %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, karl.value), y=karl.value)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/CLAMS/prolif_karl_all.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/CLAMS/prolif_karl_all.pdf", width=9.3, height=5)


# karl box - final used
patient_annotation_all <- patient_annotation_all %>% mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
groups.to.analyze <- levels(factor(patient_annotation_all$groups.to.analyze))
prolif_summary_tru <- patient_annotation_all[patient_annotation_all$clams.class == 'TRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_tru <- prolif_summary_tru %>% rename(tru = n)
prolif_summary_nontru <- patient_annotation_all[patient_annotation_all$clams.class == 'NonTRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_nontru <- prolif_summary_nontru %>% rename(nontru = n)
prolif_summary_tru_non <- merge(x=prolif_summary_nontru, y=prolif_summary_tru,
                                by="groups.to.analyze", all=TRUE)
rm(prolif_summary_nontru)
rm(prolif_summary_tru)

patient_annotation_all %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(data=patient_annotation_all, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_clams) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=groups.to.analyze, y=7.3e+05, label=prolif_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=0, y=7.3e+05, label="n =", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=0, y=2.3e+05, label="n =", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=groups.to.analyze, y=2.3e+05, label=prolif_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/CLAMS/prolif_karl_clams.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/CLAMS/prolif_karl_clams.pdf", width=9.3, height=5)

# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, karl.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, karl.value), y=karl.value, fill=clams.class)) +
#   geom_boxplot(position="identity", outlier.size = 0.3) +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Cell proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   theme(axis.text.x = element_text(angle = 60, hjust=1)) +
#   scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))

# karl box - outlier color
patient_annotation_all %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, karl.value), y=karl.value, col=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3, outlier.color = NULL) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_clams) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))

# # karl box flip
# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, karl.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, karl.value), y=karl.value, fill=clams.class)) +
#   geom_boxplot(outlier.size = 0.1) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))
# 
# # karl violin
# patient_annotation_all %>%
#   mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
#   mutate(graph.name = fct_reorder(graph.name, karl.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(graph.name, karl.value), y=karl.value, fill=clams.class)) +
#   geom_violin(position="identity") +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   theme(axis.text.x = element_text(angle = 60, hjust=1)) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))

# immune box - no CLAMS separation
patient_annotation_all %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(3e+05, 6e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all.pdf", width=9.3, height=5)

patient_annotation_all$clams.class <- NA
tru <- patient_annotation_all_p[ which(patient_annotation_all_p$clams.class == "TRU"), ] %>% pull(sample.id)
non <- patient_annotation_all_p[ which(patient_annotation_all_p$clams.class == "NonTRU"), ] %>% pull(sample.id)
patient_annotation_all <- mutate(patient_annotation_all, clams.class = 
                                   ifelse(patient_annotation_all$sample.id %in% tru, "TRU", 
                                          ifelse(patient_annotation_all$sample.id %in% non, "NonTRU", clams.class)))

# immune box - clams
patient_annotation_all %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(3e+05, 6e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams.pdf", width=9.3, height=5)

# immune box - clams - outlier color
patient_annotation_all %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value, col=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3, outlier.color = NULL) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))

#### IMMUNE WITH TUMOR CONTENT   ------

cancer.fraction.tcga <- patient_annotation_tcga[c("sample.id", "Cancer.DNA.fraction", "clams.class")]
patient_annotation_all <- left_join(patient_annotation_all, cancer.fraction.tcga)
patient_annotation_all <- patient_annotation_all %>% mutate(immune.corrected = immune.value * Cancer.DNA.fraction)

patient_annotation_all_no_nas <- patient_annotation_all[!is.na(patient_annotation_all$immune.corrected),]

# immune corrected no clams
patient_annotation_all_no_nas %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=immune.corrected)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(1e+05, 5e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_corrected.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_corrected.pdf", width=9.3, height=5)

immune_summary_tru <- patient_annotation_all_no_nas[patient_annotation_all_no_nas$clams.class == 'TRU', ] %>% 
  group_by(cancer.type) %>% tally()
immune_summary_tru <- immune_summary_tru %>% rename(tru = n)
immune_summary_nontru <- patient_annotation_all_no_nas[patient_annotation_all_no_nas$clams.class == 'NonTRU', ] %>% 
  group_by(cancer.type) %>% tally()
immune_summary_nontru <- immune_summary_nontru %>% rename(nontru = n)
immune_summary_tru_non <- merge(x=immune_summary_nontru, y=immune_summary_tru,
                                by="cancer.type", all=TRUE)
rm(immune_summary_nontru)
rm(immune_summary_tru)

patient_annotation_all_no_nas <- patient_annotation_all_no_nas %>% 
  mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
groups.to.analyze <- levels(factor(patient_annotation_all_no_nas$groups.to.analyze))

patient_annotation_all_no_nas %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=immune.corrected, fill=clams.class)) +
  geom_boxplot(position="dodge", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(breaks = c(1e+05, 5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=groups.to.analyze, y=6.5e+05, label=immune_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=-0.1, y=6.5e+05, label="n =", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=-0.1, y=0.3e+05, label="n =", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=groups.to.analyze, y=0.3e+05, label=immune_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams_corrected.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams_corrected.pdf", width=9.3, height=5)

# immune cancer.dna.fraction
patient_annotation_all_no_nas %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=Cancer.DNA.fraction)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Cancer DNA fraction", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_dna_fraction.png", width=9.3, height=5, dpi=300)

# NEW PROLIFERATION ANALYSIS PER CANCER TYPE WITH KARL VALUES --------------------------------------------------------

patient_annotation_all <- patient_annotation_all %>% mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
groups.to.analyze <- levels(factor(patient_annotation_all$groups.to.analyze))
patient_annotation_all <- mutate(patient_annotation_all, prolif.group.karl = "not classified")

for (cancer.group in groups.to.analyze) {
  print(cancer.group)
  # get only the data for that dataset + cancer type
  current.data <- subset(patient_annotation_all, groups.to.analyze == cancer.group)
  # if karl < 25th percentile, call it "Low", else "High"
  current.data$prolif.group.karl <- ifelse(current.data$karl.value < quantile(current.data$karl.value, .25), "Low", "High")
  lows <- current.data[ which(current.data$prolif.group.karl == "Low"), ] %>% pull(sample.id)
  highs <- current.data[ which(current.data$prolif.group.karl == "High"), ] %>% pull(sample.id)
  patient_annotation_all <- mutate(patient_annotation_all, prolif.group.karl = 
                                     ifelse(patient_annotation_all$sample.id %in% lows, "Low", 
                                            ifelse(patient_annotation_all$sample.id %in% highs, "High", prolif.group.karl)))
}

setwd("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/")
write.csv(patient_annotation_all, "prolif_all_samples_new.csv", row.names=FALSE)

# taking a look - boxplot
prolif_summary_high <- patient_annotation_all[patient_annotation_all$prolif.group.karl == 'High', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_low <- patient_annotation_all[patient_annotation_all$prolif.group.karl == 'Low', ] %>% 
  group_by(groups.to.analyze) %>% tally()

patient_annotation_all %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, fill=prolif.group.karl)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#26828E", "#FDE725"), name="Proliferation") +
  labs(data=patient_annotation_all, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = sign_os_prolif) +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major.x = element_blank(), axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=groups.to.analyze, y=7.3e+05, label=prolif_summary_high$n, size=3, color="#26828E", angle=90) +
  annotate("text", x=-0.1, y=7.3e+05, label="n =", size=3, color="#26828E", hjust=0) +
  annotate("text", x=-0.1, y=2.3e+05, label="n =", size=3, color="#FDE725", hjust=0) +
  annotate("text", x=groups.to.analyze, y=2.3e+05, label=prolif_summary_low$n, size=3, color="#FDE725", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/25th percentile/prolif_karl_25p.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/25th percentile/prolif_karl_25p.pdf", width=9.3, height=5)


## ANALYSIS ------

# get patient annotation tables as above

patient_annotation_gobo <- readRDS("/media/deboraholi/Data/LUND/9 THESIS/data/breast_GOBO/patient_clams.rds")
patient_annotation_scanb <- readRDS("/media/deboraholi/Data/LUND/9 THESIS/data/breast_SCAN_B/patient_clams.rds")
patient_annotation_tcga <- readRDS("/media/deboraholi/Data/LUND/9 THESIS/data/all_tcga/patient_clams.rds")

# add column with prolif.group.karl
lows <- patient_annotation_all[ which(patient_annotation_all$prolif.group.karl == "Low"), ] %>% pull(sample.id)
highs <- patient_annotation_all[ which(patient_annotation_all$prolif.group.karl == "High"), ] %>% pull(sample.id)

patient_annotation_gobo <- mutate(patient_annotation_gobo, prolif.group.karl = "not classified")
patient_annotation_gobo <- mutate(patient_annotation_gobo, prolif.group.karl = 
                                    ifelse(patient_annotation_gobo$SampleID %in% lows, "Low", 
                                           ifelse(patient_annotation_gobo$SampleID %in% highs, "High", prolif.group.karl)))

patient_annotation_scanb <- mutate(patient_annotation_scanb, prolif.group.karl = "not classified")
patient_annotation_scanb <- mutate(patient_annotation_scanb, prolif.group.karl = 
                                     ifelse(patient_annotation_scanb$rba %in% lows, "Low", 
                                            ifelse(patient_annotation_scanb$rba %in% highs, "High", prolif.group.karl)))

patient_annotation_tcga <- mutate(patient_annotation_tcga, prolif.group.karl = "not classified")
patient_annotation_tcga <- mutate(patient_annotation_tcga, prolif.group.karl = 
                                    ifelse(patient_annotation_tcga$sample_id %in% lows, "Low", 
                                           ifelse(patient_annotation_tcga$sample_id %in% highs, "High", prolif.group.karl)))

# quick chi-square
counts <- table(patient_annotation_all$clams.class, patient_annotation_all$prolif.group.karl)
counts
chisq.test(counts)
# p-value < 2.2e-16, so yes, low proliferative correlated with TRU and high proliferative with nonTRU


library(survival)
library(survminer)


## GOBO ------------

os_time_column_to_use <- "OS"
overall_survival_event_column_name <- "OSbin"

# CENSOR DATA AT SPECIFIC NUMBER OF DAYS/MONTHS/YEARS
censor_at_timepoint <- 5
number_to_cap_at <- 5
# if OS number > number_to_cap_at, set time to cap, else keep original number
patient_annotation_gobo <- mutate(patient_annotation_gobo, new_os_time_column =
                                    ifelse(
                                      get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                      as.numeric(number_to_cap_at), get(os_time_column_to_use)))
# if OS number > number_to_cap_at, set event to 0, else keep original event
patient_annotation_gobo <- mutate(patient_annotation_gobo, new_os_event_column =
                                    ifelse(
                                      get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                      0, get(overall_survival_event_column_name)))
new_os_time_column_name <- paste(os_time_column_to_use,
                                 censor_at_timepoint, sep="_")
new_os_event_column_name <- paste(overall_survival_event_column_name,
                                  censor_at_timepoint, sep="_")
names(patient_annotation_gobo)[names(patient_annotation_gobo) == "new_os_time_column"] <- new_os_time_column_name
names(patient_annotation_gobo)[names(patient_annotation_gobo) == "new_os_event_column"] <- new_os_event_column_name  

# GET THE CORRECT COLUMN FOR OS TIME (DAYS/MONTHS/YEARS, CENSORED/NOT)
time_column_to_use <- new_os_time_column_name
event_column_to_use <- new_os_event_column_name

# COX for hazard ratios - CLAMS
patient_annotation_gobo$clams_class <- factor(patient_annotation_gobo$clams_class, levels=c("TRU", "NonTRU"))
summary(coxph(formula = Surv(OS_5, OSbin_5)~clams_class, data = patient_annotation_gobo))
print(ggforest(coxph(formula = Surv(OS_5, OSbin_5)~clams_class, data = patient_annotation_gobo), main = "Hazard ratio GOBO"))

# COX for hazard ratios - prolif.group.karl
patient_annotation_gobo$prolif.group.karl <- factor(patient_annotation_gobo$prolif.group.karl, levels=c("Low", "High"))
summary(coxph(formula = Surv(OS_5, OSbin_5)~prolif.group.karl, data = patient_annotation_gobo))
print(ggforest(coxph(formula = Surv(OS_5, OSbin_5)~prolif.group.karl, data = patient_annotation_gobo), main = "Hazard ratio GOBO"))






## SCAN-B ------------

os_time_column_to_use <- "OS"
overall_survival_event_column_name <- "OSbin"

# CENSOR DATA AT SPECIFIC NUMBER OF DAYS/MONTHS/YEARS
censor_at_timepoint <- 5
number_to_cap_at <- 5
# if OS number > number_to_cap_at, set time to cap, else keep original number
patient_annotation_scanb <- mutate(patient_annotation_scanb, new_os_time_column =
                                     ifelse(
                                       get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                       as.numeric(number_to_cap_at), get(os_time_column_to_use)))
# if OS number > number_to_cap_at, set event to 0, else keep original event
patient_annotation_scanb <- mutate(patient_annotation_scanb, new_os_event_column =
                                     ifelse(
                                       get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                       0, get(overall_survival_event_column_name)))
new_os_time_column_name <- paste(os_time_column_to_use,
                                 censor_at_timepoint, sep="_")
new_os_event_column_name <- paste(overall_survival_event_column_name,
                                  censor_at_timepoint, sep="_")
names(patient_annotation_scanb)[names(patient_annotation_scanb) == "new_os_time_column"] <- new_os_time_column_name
names(patient_annotation_scanb)[names(patient_annotation_scanb) == "new_os_event_column"] <- new_os_event_column_name  

# GET THE CORRECT COLUMN FOR OS TIME (DAYS/MONTHS/YEARS, CENSORED/NOT)
time_column_to_use <- new_os_time_column_name
event_column_to_use <- new_os_event_column_name

# COX for hazard ratios - CLAMS
patient_annotation_scanb$clams_class <- factor(patient_annotation_scanb$clams_class, levels=c("TRU", "NonTRU"))
summary(coxph(formula = Surv(OS_5, OSbin_5)~clams_class, data = patient_annotation_scanb))
print(ggforest(coxph(formula = Surv(OS_5, OSbin_5)~clams_class, data = patient_annotation_scanb), main = "Hazard ratio SCAN-B"))

# COX for hazard ratios - prolif.group.karl
patient_annotation_scanb$prolif.group.karl <- factor(patient_annotation_scanb$prolif.group.karl, levels=c("Low", "High"))
summary(coxph(formula = Surv(OS_5, OSbin_5)~prolif.group.karl, data = patient_annotation_scanb))
print(ggforest(coxph(formula = Surv(OS_5, OSbin_5)~prolif.group.karl, data = patient_annotation_scanb), main = "Hazard ratio SCAN-B"))





## TCGA ------------

os_time_column_to_use <- "OS.time"
overall_survival_time_column_name <- "OS.time"
overall_survival_event_column_name <- "OS"
overall_survival_time_column_unit_table <- "days"
overall_survival_time_column_unit_desired <- "years"

# MUTATE patient_annotation_tcga SO OS TIME IS IN MONTHS OR YEARS INSTEAD OF DAYS IF NEEDED
if (overall_survival_time_column_unit_table == overall_survival_time_column_unit_desired) {
  os_time_column_to_use <- overall_survival_time_column_name
} else if (overall_survival_time_column_unit_table == 'days' & 
           overall_survival_time_column_unit_desired == 'years') {
  patient_annotation_tcga <- mutate(patient_annotation_tcga, os_time_years = get(overall_survival_time_column_name) / 365.2425)
  new_os_time_column_name <- paste0(overall_survival_time_column_name, "_years")
  names(patient_annotation_tcga)[names(patient_annotation_tcga) == "os_time_years"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else if (overall_survival_time_column_unit_table == 'months' & 
           overall_survival_time_column_unit_desired == 'years') {
  patient_annotation_tcga <- mutate(patient_annotation_tcga, os_time_years = get(overall_survival_time_column_name) / 12)
  new_os_time_column_name <- paste0(overall_survival_time_column_name, "_years")
  names(patient_annotation_tcga)[names(patient_annotation_tcga) == "os_time_years"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else if (overall_survival_time_column_unit_table == 'days' & 
           overall_survival_time_column_unit_desired == 'months') {
  patient_annotation_tcga <- mutate(patient_annotation_tcga, os_time_months = get(overall_survival_time_column_name) / 30.436875)
  new_os_time_column_name <- paste0(overall_survival_time_column_name, "_months")
  names(patient_annotation_tcga)[names(patient_annotation_tcga) == "os_time_months"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else {
  print('Something went wrong when changing the time unit')
}


# CENSOR DATA AT SPECIFIC NUMBER OF DAYS/MONTHS/YEARS
censor_at_timepoint <- 5
number_to_cap_at <- 5
# if OS number > number_to_cap_at, set time to cap, else keep original number
patient_annotation_tcga <- mutate(patient_annotation_tcga, new_os_time_column =
                                    ifelse(
                                      get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                      as.numeric(number_to_cap_at), get(os_time_column_to_use)))
# if OS number > number_to_cap_at, set event to 0, else keep original event
patient_annotation_tcga <- mutate(patient_annotation_tcga, new_os_event_column =
                                    ifelse(
                                      get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                      0, get(overall_survival_event_column_name)))
new_os_time_column_name <- paste(os_time_column_to_use,
                                 censor_at_timepoint, sep="_")
new_os_event_column_name <- paste(overall_survival_event_column_name,
                                  censor_at_timepoint, sep="_")
names(patient_annotation_tcga)[names(patient_annotation_tcga) == "new_os_time_column"] <- new_os_time_column_name
names(patient_annotation_tcga)[names(patient_annotation_tcga) == "new_os_event_column"] <- new_os_event_column_name  

# GET THE CORRECT COLUMN FOR OS TIME (DAYS/MONTHS/YEARS, CENSORED/NOT)
time_column_to_use <- new_os_time_column_name
event_column_to_use <- new_os_event_column_name

# COX for hazard ratios - CLAMS
patient_annotation_tcga$clams_class <- factor(patient_annotation_tcga$clams_class, levels=c("TRU", "NonTRU"))
summary(coxph(formula = Surv(get(time_column_to_use), get(event_column_to_use))~clams_class, data = patient_annotation_tcga))
print(ggforest(coxph(formula = Surv(get(time_column_to_use), get(event_column_to_use))~clams_class, data = patient_annotation_tcga), main = "Hazard ratio TCGA BRCA"))

# COX for hazard ratios - prolif.group.karl
patient_annotation_tcga$prolif.group.karl <- factor(patient_annotation_tcga$prolif.group.karl, levels=c("Low", "High"))
groups.to.analyze <- levels(factor(patient_annotation_tcga$cancer.type))

for (type_to_analyze in groups.to.analyze) {
  print(type_to_analyze)
  type_subset <- subset(patient_annotation_tcga, cancer.type == type_to_analyze)
  type_subset$prolif.group.karl <- factor(type_subset$prolif.group.karl, levels=c("Low", "High"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column_to_use]], 
                              event=type_subset[[event_column_to_use]])
  print(summary(coxph(formula = overall_surv_object~prolif.group.karl, data = type_subset)))
  print(ggforest(coxph(formula = overall_surv_object~prolif.group.karl, data = type_subset), main = paste("Hazard ratio", type_to_analyze)))
}



# OVERALL SURVIVAL ----

# FOR ONLY ONE CANCER TYPE IN THE DATASET - GOBO OR SCAN-B

patient_annotation <- patient_annotation_gobo
overall_survival_name <- "BRCA - GOBO"


patient_annotation <- patient_annotation_scanb
overall_survival_name <- "BRCA - SCAN-B"

# analysis
patient_annotation$prolif.group.karl <- factor(patient_annotation$prolif.group.karl, levels=c("Low", "High"))
# make the object
overall_surv_object <- Surv(time=patient_annotation[[time_column_to_use]], 
                            event=patient_annotation[[event_column_to_use]])
# separating estimates by Low and High proliferative
fit_prolif <- survfit(overall_surv_object~prolif.group.karl, data=patient_annotation)
# testing by clams class
print(survdiff(Surv(patient_annotation[[time_column_to_use]], 
                    patient_annotation[[event_column_to_use]])~patient_annotation$prolif.group.karl,))

# plot
# simple for comparing info and making sure names are correct
print(ggsurvplot(fit_prolif, data=patient_annotation, title="Overall Survival", risk.table=TRUE))
# good looking one
current_plot <- ggsurvplot(fit_prolif, data=patient_annotation,
                           palette=c("#26828E", "#FDE725"),
                           title=paste0("Overall Survival (", overall_survival_name, ")"), 
                           xlab=paste0("Time (years)"),
                           censor.shape=124, censor.size=3,
                           pval=TRUE, pval.coord=c(0,0.1),
                           surv.median.line="hv",
                           risk.table=TRUE,
                           risk.table.fontsize = 4,
                           tables.theme = theme_survminer(font.main = 14),
                           legend="none", legend.title="Proliferation",
                           legend.labs=c("Low", "High"))
print(current_plot)
current_filename <- "OS 5y GOBO.png" # GOBO
current_filename <- "OS 5y SCANB.png" # SCAN-B
ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)


# FOR TCGA

patient_annotation <- patient_annotation_tcga
patient_annotation$prolif.group.karl <- factor(patient_annotation$prolif.group.karl, levels=c("Low", "High"))
groups.to.analyze <- levels(factor(patient_annotation$cancer.type))

groups.to.analyze <- "LUAD"

setwd("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/25th percentile/OS 5y/")

for (type_to_analyze in groups.to.analyze) {
  print(type_to_analyze)
  type_subset <- subset(patient_annotation, cancer.type == type_to_analyze)
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column_to_use]], 
                              event=type_subset[[event_column_to_use]])
  # separating estimates by TRU or NonTRU
  fit_prolif <- survfit(overall_surv_object~prolif.group.karl, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column_to_use]], 
                      type_subset[[event_column_to_use]])~type_subset$prolif.group.karl,))
  # plot
  # simple for comparing info and making sure names are correct
  #print(ggsurvplot(fit_prolif, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  current_plot <- ggsurvplot(fit_prolif, data=type_subset,
                             palette=c("#26828E", "#FDE725"),
                             title=paste0("Overall Survival (", type_to_analyze, ")"), 
                             xlab=paste0("Time (years)"),
                             censor.shape=124, censor.size=3,
                             pval=TRUE, pval.coord=c(0,0.1),
                             surv.median.line="hv",
                             risk.table=TRUE,
                             risk.table.fontsize = 4,
                             tables.theme = theme_survminer(font.main = 14),
                             legend="none", legend.title="Proliferation",
                             legend.labs=c("Low", "High"))
  print(current_plot)
  current_filename <- paste0("OS 5y ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}


# ALL PROLIFERATION PLOTS TOGETHER ------

# make a summary of general karl.value
patient_annotation_prolif_all <- patient_annotation_all %>% 
  group_by(groups.to.analyze) %>% 
  summarise_at("karl.value", median)

groups.to.analyze <- levels(factor(patient_annotation_all$groups.to.analyze))
patient_annotation_prolif_all$groups.to.analyze <- factor(patient_annotation_prolif_all$groups.to.analyze)

patient_annotation_prolif_all %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value)) +
  geom_point(color = "gray", shape=15, size=3) +
  geom_line(aes(x=groups.to.analyze, y=karl.value, group=1), color = "gray") +
  # karl.value by low/high
  geom_point(data= patient_annotation_all_new %>% 
               group_by(groups.to.analyze, prolif.group.karl) %>% 
               summarise_at("karl.value", median),
             aes(x=groups.to.analyze, y=karl.value, group=prolif.group.karl, color=prolif.group.karl),
             size=2) +
  geom_line(data=patient_annotation_all_new %>% 
              group_by(groups.to.analyze, prolif.group.karl) %>% 
              summarise_at("karl.value", median),
            aes(x=groups.to.analyze, y=karl.value, group=prolif.group.karl, color=prolif.group.karl)) +
  # karl.value by clams
  geom_point(data=patient_annotation_all %>% 
               group_by(groups.to.analyze, clams.class) %>% 
               summarise_at("karl.value", median),
             aes(x=groups.to.analyze, y=karl.value, group=clams.class, color=clams.class),
             shape=18, size=3) +
  geom_line(data=patient_annotation_all %>% 
              group_by(groups.to.analyze, clams.class) %>% 
              summarise_at("karl.value", median),
            aes(x=groups.to.analyze, y=karl.value, group=clams.class, color=clams.class)) +
  scale_color_manual(values=c("#6DCD59", "#FDE725", "deepskyblue3", "darkorange1")) +
  theme_minimal() +
  labs(y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        legend.position=c(0.15,0.75)) +
  scale_y_continuous(breaks = c(3.5e+05, 6e+05), labels = c('Lower', 'Higher'))


ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/prolif_karl_lines.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Proliferation Karlsson/prolif_karl_lines.pdf", width=9.3, height=5)

## FURTHER CHARACTERIZE KIDNEY SAMPLES -----------

# read in the table with kidney classifications
kidney_class <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/data/kidney_clams_taxonomy.csv")
kidney_class <- kidney_class %>% mutate(cancer.type = substr(sample_id, 1, 4))
kidney_class$cancer.type <- factor(kidney_class$cancer.type)
kidney_class$clams.class <- factor(kidney_class$clams.class)
kidney_class$taxonomy_published <- factor(kidney_class$taxonomy_published)
kidney_taxonomy <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/data/kidney_taxonomy.csv")
kidney_class <- left_join(kidney_class, kidney_taxonomy, by=c("taxonomy_published"="predom_enriched"))

kidney_class_summary <- kidney_class %>% count(cancer.type, clams.class, predom, taxonomy_published)

addmargins(table(kidney_class$clams.class, kidney_class$predom, kidney_class$cancer.type, useNA = "ifany"))
addmargins(table(kidney_class$clams.class, kidney_class$taxonomy_published, kidney_class$cancer.type, useNA = "ifany"))

for (cancer.type in levels(kidney_class$cancer.type)) {
  type.subset <- kidney_class[kidney_class$cancer.type == cancer.type, ]
  print(ggplot(type.subset, aes(x=clams.class, fill=taxonomy_published)) +
          geom_bar(position = 'stack') +
          theme_minimal() +
          ggtitle(cancer.type))
}


## TCGA MEDIAN OS TIMES -----------

tcga_median_OS <- patient_annotation_tcga %>% group_by(cancer.type) %>% summarise_at(vars(OS.time), median, na.rm = TRUE)
tcga_median_OS <- tcga_median_OS %>% rename(OS.time.days = OS.time)
tcga_median_OS <- tcga_median_OS %>% mutate(OS.time.months = round(OS.time.days * 0.0328549112, digits=1))
tcga_median_OS <- tcga_median_OS %>% mutate(OS.time.years = round(OS.time.days * 0.00273790926, digits=2))

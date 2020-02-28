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




# PLOTS  ------------------------------------------------------------------------

# if starting from here
# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# CLAMS RESULTS
# IMMUNE ANALYSIS RESULTS
# PLOTS NAMES

patient_annotation_immune_clams <- left_join(patient_annotation_immune, patient_annotation_clams, 
                                             by=c("sample.id", "cancer.type", "dataset"))


# COMPARE MAX AND SUM AS WAY OF REMOVING DUPLICATES  ---------------------------

# GOBO

gobo_subset <- patient_annotation_immune_clams %>% subset(dataset == "GOBO")

# Trying to see with points
gobo_subset %>%
  ggplot(aes(x=immune.value.max, y=immune.value.sum, color=clams.class)) + 
  geom_point() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1")) +
  theme_minimal() +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"))
# Correlation test
cor.test(gobo_subset$immune.value.max, gobo_subset$immune.value.sum) #0.9940132

# using max is fine

ggplot(gobo_subset %>% subset(clams.class == "TRU"), aes(x=immune.value.max, y=immune.value.sum, color=clams.class)) + 
  geom_point() +
  scale_color_manual(values=c("darkorange1")) +
  theme_minimal() +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"))
ggplot(gobo_subset %>% subset(clams.class == "NonTRU"), aes(x=immune.value.max, y=immune.value.sum, color=clams.class)) + 
  geom_point() +
  scale_color_manual(values=c("deepskyblue3")) +
  theme_minimal() +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"))


# PLOTS  ---------------------------------------------------

# rename desired column to be used in the plots (sum or max)
patient_annotation_immune_clams <- patient_annotation_immune_clams %>% 
                                      rename(immune.value = immune.value.max)

# boxplot, no CLAMS separation
patient_annotation_immune_clams %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all.pdf", width=9.3, height=5)


# boxplot by CLAMS
patient_annotation_immune_clams %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(3e+05, 6e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams.pdf", width=9.3, height=5)

# boxplot by CLAMS wih colored outliers
patient_annotation_immune_clams %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.value, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.value), y=immune.value, col=clams.class, alpha = 0.2)) +
  geom_boxplot(position="identity", outlier.size = 0.3, outlier.color = NULL) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_alpha(guide = 'none') +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher'))



#### IMMUNE CORRECTED BY TUMOR CONTENT   ---------------------------------------
# only possible for the TCGA dataset


# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# CLAMS RESULTS
# IMMUNE ANALYSIS RESULTS
# (1) PATIENT ANNOTATION - only for TCGA

# if needed
patient_annotation_immune_clams <- left_join(patient_annotation_immune, patient_annotation_clams, 
                                             by=c("sample.id", "cancer.type", "dataset"))


cancer_fraction <- patient_annotation_tcga[c("sample_id", "Cancer.DNA.fraction")]
patient_annotation_fraction <- inner_join(patient_annotation_immune_clams, cancer_fraction,
                                                       by = c("sample.id" = "sample_id"))
patient_annotation_fraction <- patient_annotation_fraction[!is.na(patient_annotation_fraction$Cancer.DNA.fraction),]
patient_annotation_fraction <- patient_annotation_fraction %>% 
                                      mutate(immune.corrected = immune.value * Cancer.DNA.fraction)

# boxplot, no CLAMS separation
patient_annotation_fraction %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=immune.corrected)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(1e+05, 5e+05), labels = c('Lower', 'Higher'))

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_corrected.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_all_corrected.pdf", width=9.3, height=5)


# boxplot by CLAMS with sample sizes
# get sample sizes
immune_summary_tru <- patient_annotation_fraction[patient_annotation_fraction$clams.class == 'TRU', ] %>% 
  group_by(cancer.type) %>% tally()
immune_summary_tru <- immune_summary_tru %>% rename(tru = n)
immune_summary_nontru <- patient_annotation_fraction[patient_annotation_fraction$clams.class == 'NonTRU', ] %>% 
  group_by(cancer.type) %>% tally()
immune_summary_nontru <- immune_summary_nontru %>% rename(nontru = n)
immune_summary_tru_non <- merge(x=immune_summary_nontru, y=immune_summary_tru,
                                by="cancer.type", all=TRUE)
rm(immune_summary_nontru)
rm(immune_summary_tru)

patient_annotation_fraction <- patient_annotation_fraction %>% 
  mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
groups.to.analyze <- levels(factor(patient_annotation_fraction$groups.to.analyze))

# fill coloring
patient_annotation_fraction %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=immune.corrected, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(1e+05, 5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=groups.to.analyze, y=6.5e+05, label=immune_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=-0.1, y=6.5e+05, label="n =", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=-0.1, y=0.3e+05, label="n =", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=groups.to.analyze, y=0.3e+05, label=immune_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

# outline coloring
patient_annotation_fraction %>%
  mutate(graph.name = paste(cancer.type, dataset, sep = '_')) %>%
  mutate(graph.name = fct_reorder(graph.name, immune.corrected, .fun='median')) %>%
  ggplot(aes(x=reorder(graph.name, immune.corrected), y=immune.corrected, col=clams.class, alpha = 0.2)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(y = "Immune signature", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_y_continuous(breaks = c(1e+05, 5e+05), labels = c('Lower', 'Higher')) +
  scale_alpha(guide = 'none') +
  annotate("text", x=groups.to.analyze, y=6.5e+05, label=immune_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=-0.1, y=6.5e+05, label="n =", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=-0.1, y=0.3e+05, label="n =", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=groups.to.analyze, y=0.3e+05, label=immune_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams_corrected.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/2_immune/immune_clams_corrected.pdf", width=9.3, height=5)


# any stats on Cancer.DNA.fraction by cancer.type?

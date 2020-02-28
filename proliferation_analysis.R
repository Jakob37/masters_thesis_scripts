#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# SCRIPT TO PREPARE PAN CANCER DATA, GET PROLIFERATION MEASURE AND PLOT
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

setwd("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/")



# GET GENE SYMBOLS AND LIST  -------------------------------------

# GOBO
reporter_to_symbol_gobo <- gex_matrix_gobo[,1:2]
reporter_to_symbol_gobo <- rownames_to_column(data.frame(reporter_to_symbol_gobo), var = "reporterId")
reporter_to_symbol_gobo <- left_join(reporter_to_symbol_gobo, gene_table_gobo[,c("reporterId","geneSymbol")])
gene_symbol_column_gobo <- reporter_to_symbol_gobo["geneSymbol"]

rm(gene_table_gobo)

length(gene_symbol_column_gobo$geneSymbol) #18223
length(unique(gene_symbol_column_gobo$geneSymbol)) #12201

genes_gobo <- unique(gene_symbol_column_gobo$geneSymbol)

# SCAN-B
# row.names already are the gene symbols, no duplicate genes
reporter_to_symbol_scanb <- gex_matrix_scanb[,1:2]
reporter_to_symbol_scanb <- rownames_to_column(data.frame(gex_matrix_scanb), var = "geneSymbol")
gene_symbol_column_scanb <- reporter_to_symbol_scanb["geneSymbol"]

length(gene_symbol_column_scanb$geneSymbol) #19102
length(unique(gene_symbol_column_scanb$geneSymbol)) #19102

genes_scanb <- unique(gene_symbol_column_scanb$geneSymbol)


#TCGA
reporter_to_symbol_tcga <- gex_matrix_tcga[,1:2]
reporter_to_symbol_tcga <- rownames_to_column(data.frame(reporter_to_symbol_tcga), var = "reporterId")
reporter_to_symbol_tcga <- left_join(reporter_to_symbol_tcga, gene_table_tcga[,c("ENSG", "SYMBOL")], by=c("reporterId" = "ENSG"))
reporter_to_symbol_tcga <- reporter_to_symbol_tcga %>% rename(geneSymbol=SYMBOL)
gene_symbol_column_tcga <- reporter_to_symbol_tcga["geneSymbol"]

rm(gene_table_tcga)

length(gene_symbol_column_tcga$geneSymbol) #19676
length(unique(gene_symbol_column_tcga$geneSymbol)) #19668

genes_tcga <- unique(gene_symbol_column_tcga$geneSymbol)



# GET COMMON GENE SYMBOLS  ---------------------------------

common_genes <- intersect(intersect(genes_gobo, genes_scanb), genes_tcga) 

length(common_genes) #10,916
length(unique(common_genes)) #10,916

rm(genes_gobo)
rm(genes_scanb)
rm(genes_tcga)



# INPUT MODULES ------------------------------------------

# Fredlund
mod_prolif_fredlund <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/Fredlund.csv")
mod_prolif_fredlund <- mod_prolif_fredlund$geneSymbol
length(mod_prolif_fredlund) #20
mod_prolif_fredlund

mod_prolif_fredlund_in_common <- mod_prolif_fredlund[mod_prolif_fredlund %in% common_genes]
length(mod_prolif_fredlund_in_common) #19
mod_prolif_fredlund_in_common

rm(mod_prolif_fredlund)

# Result
write.table(mod_prolif_fredlund_in_common, 
            file='/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/prolif_fred_genes_in_common.txt', 
            quote=FALSE, row.names=FALSE, col.names=FALSE)


# Karlsson
mod_prolif_karlsson <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/Karlsson.csv")
mod_prolif_karlsson <- mod_prolif_karlsson$geneSymbol
length(mod_prolif_karlsson) #93
mod_prolif_karlsson

mod_prolif_karlsson_in_common <- mod_prolif_karlsson[mod_prolif_karlsson %in% common_genes]
length(mod_prolif_karlsson_in_common) #75
mod_prolif_karlsson_in_common

rm(mod_prolif_karlsson)

# Result
write.table(mod_prolif_karlsson_in_common, 
            file='/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/modules/prolif_karl_genes_in_common.txt', 
            quote=FALSE, row.names=FALSE, col.names=FALSE)



# MODULE ANALYSIS  ------------------------------------

# GET SIGNATURE FUNCTION FROM FILE
source("/media/deboraholi/Data/LUND/9 THESIS/src/functions.R")


# GOBO

# Fredlund
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_gobo))) {
  sample_column <- gex_matrix_gobo[,i,drop=FALSE]
  prolif_value_fred <- get_signature_value(gene_symbol_column_gobo, sample_column, mod_prolif_fredlund_in_common, max)
  patient_ids_all_datasets$fred.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_fred
  rm(sample_column)
  if (i %% 100 == 0) { print(i) }
}

# Karlsson
for (i in 1:length(colnames(gex_matrix_gobo))) {
  sample_column <- gex_matrix_gobo[,i,drop=FALSE]
  prolif_value_karl <- get_signature_value(gene_symbol_column_gobo, sample_column, mod_prolif_karlsson_in_common, max)
  patient_ids_all_datasets$karl.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_karl
  rm(sample_column)
  if (i %% 100 == 0) { print(i) }
}

rm(gex_matrix_gobo)


# SCAN-B

# Fredlund
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_scanb))) {
  sample_column <- gex_matrix_scanb[,i,drop=FALSE]
  prolif_value_fred <- get_signature_value(gene_symbol_column_scanb, sample_column, mod_prolif_fredlund_in_common, max)
  patient_ids_all_datasets$fred.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_fred
  if (i %% 100 == 0) { print(i) }
}

# Karlsson
for (i in 1:length(colnames(gex_matrix_scanb))) {
  sample_column <- gex_matrix_scanb[,i,drop=FALSE]
  prolif_value_karl <- get_signature_value(gene_symbol_column_scanb, sample_column, mod_prolif_karlsson_in_common, max)
  patient_ids_all_datasets$karl.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_karl
  if (i %% 100 == 0) { print(i) }
}

rm(gex_matrix_scanb)


# TCGA

# Fredlund
# with max as the way of removing duplicate genes (only important for GOBO and TCGA tables)
for (i in 1:length(colnames(gex_matrix_tcga))) {
  sample_column <- gex_matrix_tcga[,i,drop=FALSE]
  prolif_value_fred <- get_signature_value(gene_symbol_column_tcga, sample_column, mod_prolif_fredlund_in_common, max)
  patient_ids_all_datasets$fred.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_fred
  if (i %% 200 == 0) { print(i) }
}

# Karlsson
for (i in 1:length(colnames(gex_matrix_tcga))) {
  sample_column <- gex_matrix_tcga[,i,drop=FALSE]
  prolif_value_karl <- get_signature_value(gene_symbol_column_tcga, sample_column, mod_prolif_karlsson_in_common, max)
  patient_ids_all_datasets$karl.value[patient_ids_all_datasets$sample.id == colnames(sample_column)[1]] <- prolif_value_karl
  if (i %% 200 == 0) { print(i) }
}

rm(gex_matrix_tcga)


# RESULT
# Save final table
write.csv(patient_ids_all_datasets, 
          "/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/prolif_all_samples.csv", row.names=FALSE)



# PLOTS  --------------------------------------------------------

# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# CLAMS RESULTS
# PROLIFERATION ANALYSIS RESULTS
# PLOTS NAMES

patient_annotation_prolif_clams <- left_join(patient_annotation_prolif, patient_annotation_clams, 
                                             by=c("sample.id", "cancer.type", "dataset"))
patient_annotation_prolif_clams <- patient_annotation_prolif_clams %>% 
                                    mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
#patient_annotation_prolif_clams$groups.to.analyze <- factor(patient_annotation_prolif_clams$groups.to.analyze)

# # Fredlund
# # boxplot, CLAMS classes
# patient_annotation_prolif_clams %>%
#   mutate(groups.to.analyze = fct_reorder(groups.to.analyze, fred.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(groups.to.analyze, fred.value), y=fred.value, fill=clams.class)) +
#   geom_boxplot(position="identity", outlier.size = 0.3) +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(y = "Proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
#         axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.ticks.x.bottom=element_line(color="grey50")) +
#   scale_y_continuous(breaks = c(50000, 150000), labels = c('Lower', 'Higher'))

# Karlsson
# boxplot, no CLAMS separation
patient_annotation_prolif_clams %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value)) +
  geom_boxplot(position="identity", outlier.size = 0.3, fill = "gray") +
  theme_minimal() +
  labs(y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"), 
        plot.margin = margin(t=10,r=15,b=10,l=10)) +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  coord_cartesian(clip = 'off') +
  geom_text(data = patient_annotation_prolif_clams %>% group_by(groups.to.analyze) %>% tally(),
            aes(x=groups.to.analyze, y=7.3e+05, label=n), size=3, color="grey50", angle=90) +
  annotate("text", x=34.5, y=7.3e+05, label="= n", size=3, color="grey50", hjust=0)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/prolif_karl_all.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/prolif_karl_all.pdf", width=9.3, height=5)


# boxplot, with CLAMS classification and sample sizes
# create table to annotate sample sizes
prolif_summary_tru <- patient_annotation_prolif_clams[patient_annotation_prolif_clams$clams.class == 'TRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_tru <- prolif_summary_tru %>% rename(tru = n)
prolif_summary_nontru <- patient_annotation_prolif_clams[patient_annotation_prolif_clams$clams.class == 'NonTRU', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_nontru <- prolif_summary_nontru %>% rename(nontru = n)
prolif_summary_tru_non <- merge(x=prolif_summary_nontru, y=prolif_summary_tru,
                                by="groups.to.analyze", all=TRUE)
rm(prolif_summary_nontru)
rm(prolif_summary_tru)

patient_annotation_prolif_clams %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, fill=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(data=patient_annotation_prolif_clams, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_clams) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=10,r=20,b=10,l=10)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=7.3e+05, label=prolif_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=34.5, y=7.32e+05, label="= n", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=29.5, y=2.22e+05, label="= n", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=2.2e+05, label=prolif_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/CLAMS/prolif_karl_clams.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/CLAMS/prolif_karl_clams.pdf", width=9.3, height=5)
  

# boxplot, outlier color
patient_annotation_prolif_clams %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, col=clams.class)) +
  geom_boxplot(position="identity", outlier.size = 0.3, outlier.color = NULL) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
  labs(data=patient_annotation_prolif_clams, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_clams) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=10,r=20,b=10,l=10)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=7.3e+05, label=prolif_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
  annotate("text", x=34.5, y=7.32e+05, label="= n", size=3, color="deepskyblue3", hjust=0) +
  annotate("text", x=29.5, y=2.22e+05, label="= n", size=3, color="darkorange1", hjust=0) +
  annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=2.2e+05, label=prolif_summary_tru_non$tru, size=3, color="darkorange1", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/CLAMS/prolif_karl_clams_outlier.pdf", width=9.3, height=5)



# # violin plot, CLAMS classification, sample sizes
# patient_annotation_prolif_clams %>%
#   mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
#   ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, fill=clams.class)) +
#   geom_violin(position="identity") +
#   theme_minimal() +
#   scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS") +
#   labs(data=patient_annotation_prolif_clams, y = "Cell proliferation", x = NULL) +
#   scale_x_discrete(labels = acronym_ast_os_clams) +
#   theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
#         axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.ticks.x.bottom=element_line(color="grey50"),
#         plot.margin = margin(t=10,r=20,b=10,l=10)) +
#   coord_cartesian(clip = 'off') +
#   scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
#   annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=7.3e+05, label=prolif_summary_tru_non$nontru, size=3, color="deepskyblue3", angle=90) +
#   annotate("text", x=34.5, y=7.32e+05, label="= n", size=3, color="deepskyblue3", hjust=0) +
#   annotate("text", x=29.5, y=2.22e+05, label="= n", size=3, color="darkorange1", hjust=0) +
#   annotate("text", x=prolif_summary_tru_non$groups.to.analyze, y=2.2e+05, label=prolif_summary_tru_non$tru, size=3, color="darkorange1", angle=90)





# NEW PROLIFERATION ANALYSIS PER CANCER TYPE WITH KARL VALUES --------------------------------------------------------

patient_annotation_prolif <- patient_annotation_prolif %>% mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))
groups.to.analyze <- levels(factor(patient_annotation_prolif$groups.to.analyze))
patient_annotation_prolif$prolif.group.karl <- "NA"

for (cancer.group in groups.to.analyze) {
  print(cancer.group)
  # get only the data for that dataset + cancer type
  current.data <- subset(patient_annotation_prolif, groups.to.analyze == cancer.group)
  # if karl.value < 25th percentile, call it "Low", else "High"
  current.data$prolif.group.karl <- ifelse(current.data$karl.value < quantile(current.data$karl.value, .25), "Low", "High")
  lows <- current.data[ which(current.data$prolif.group.karl == "Low"), ] %>% pull(sample.id)
  highs <- current.data[ which(current.data$prolif.group.karl == "High"), ] %>% pull(sample.id)
  patient_annotation_prolif$prolif.group.karl <- ifelse(patient_annotation_prolif$sample.id %in% lows, "Low",
                                                  ifelse(patient_annotation_prolif$sample.id %in% highs, "High", 
                                                         patient_annotation_prolif$prolif.group.karl))
}

# RESULT
write.csv(patient_annotation_prolif, 
          "/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/prolif_all_samples.csv", row.names=FALSE)

# LOAD RESULTS FROM data_input.R  -----------------------------------------------

# PROLIFERATION ANALYSIS RESULTS
# PLOTS NAMES

# taking a look - boxplot
prolif_summary_high <- patient_annotation_prolif[patient_annotation_prolif$prolif.group.karl == 'High', ] %>% 
  group_by(groups.to.analyze) %>% tally()
prolif_summary_low <- patient_annotation_prolif[patient_annotation_prolif$prolif.group.karl == 'Low', ] %>% 
  group_by(groups.to.analyze) %>% tally()

patient_annotation_prolif %>%
  mutate(groups.to.analyze = fct_reorder(groups.to.analyze, karl.value, .fun='median')) %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value, fill=prolif.group.karl)) +
  geom_boxplot(position="identity", outlier.size = 0.3) +
  theme_minimal() +
  scale_fill_manual(values=c("#26828E", "#FDE725"), name="Proliferation") +
  labs(data=patient_annotation_prolif, y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_ast_os_prolif_highlow) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        plot.margin = margin(t=0,r=15,b=0,l=0)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(breaks = c(3e+05, 6.5e+05), labels = c('Lower', 'Higher')) +
  annotate("text", x=groups.to.analyze, y=7.3e+05, label=prolif_summary_high$n, size=3, color="#26828E", angle=90) +
  annotate("text", x=34.5, y=7.3e+05, label="= n", size=3, color="#26828E", hjust=0) +
  annotate("text", x=34.5, y=2.3e+05, label="= n", size=3, color="#FDE725", hjust=0) +
  annotate("text", x=groups.to.analyze, y=2.3e+05, label=prolif_summary_low$n, size=3, color="#FDE725", angle=90)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/25th percentile/prolif_karl_25p.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/25th percentile/prolif_karl_25p.pdf", width=9.3, height=5)



# ALL PROLIFERATION PLOTS TOGETHER -------------------------------------------

# make a summary of general karl.value
patient_annotation_prolif_karl_median <- patient_annotation_prolif_clams %>% 
                                  group_by(groups.to.analyze) %>% 
                                    summarise_at("karl.value", median)

groups.to.analyze <- levels(factor(patient_annotation_prolif_clams$groups.to.analyze))
patient_annotation_prolif_karl_median$groups.to.analyze <- factor(patient_annotation_prolif_karl_median$groups.to.analyze)

patient_annotation_prolif_karl_median %>%
  ggplot(aes(x=reorder(groups.to.analyze, karl.value), y=karl.value)) +
  geom_point(color = "gray", shape=15, size=3) +
  geom_line(aes(x=groups.to.analyze, y=karl.value, group=1), color = "gray") +
  # karl.value by low/high
  geom_point(data= patient_annotation_prolif_clams %>% 
               group_by(groups.to.analyze, prolif.group.karl) %>% 
               summarise_at("karl.value", median),
             aes(x=groups.to.analyze, y=karl.value, group=prolif.group.karl, color=prolif.group.karl),
             size=2) +
  geom_line(data=patient_annotation_prolif_clams %>% 
              group_by(groups.to.analyze, prolif.group.karl) %>% 
              summarise_at("karl.value", median),
            aes(x=groups.to.analyze, y=karl.value, group=prolif.group.karl, color=prolif.group.karl)) +
  # karl.value by clams
  geom_point(data=patient_annotation_prolif_clams %>% 
               group_by(groups.to.analyze, clams.class) %>% 
               summarise_at("karl.value", median),
             aes(x=groups.to.analyze, y=karl.value, group=clams.class, color=clams.class),
             shape=18, size=3) +
  geom_line(data=patient_annotation_prolif_clams %>% 
              group_by(groups.to.analyze, clams.class) %>% 
              summarise_at("karl.value", median),
            aes(x=groups.to.analyze, y=karl.value, group=clams.class, color=clams.class)) +
  scale_color_manual(values=c("#6DCD59", "#FDE725", "deepskyblue3", "darkorange1"), name="Group",
                     labels=c("75% most proliferative", "25% less proliferative", "CLAMS NonTRU", "CLAMS TRU")) +
  theme_minimal() +
  labs(y = "Cell proliferation", x = NULL) +
  scale_x_discrete(labels = acronym_only) +
  theme(axis.line.x = element_line(color="grey50"), axis.line.y = element_line(color="grey50"),
        axis.text.x = element_text(angle = 60, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50"),
        legend.position=c(0.15,0.75)) +
  scale_y_continuous(breaks = c(3.5e+05, 6e+05), labels = c('Lower', 'Higher'))
  

ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/prolif_karl_lines.png", width=9.3, height=5, dpi=300)
ggsave("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/prolif_karl_lines.pdf", width=9.3, height=5)

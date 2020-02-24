############################
#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO GET SUMMARY DATA AND CREATE CROSS-DATA PLOTS
#
############################

# CLEAR ENVIRONMENT ---------------------------------------------------------------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ---------------------------------------------------------------------------------------------------

library(tidyverse)


# INPUT DATA ---------------------------------------------------------------------------------------------------

data_directory_path <- "/media/deboraholi/Data/LUND/9 THESIS/0_clams/"

# path to data_for_summary.csv
datasets_to_use <- c("all_tcga/data_for_summary.csv", 
                  "breast_GOBO/data_for_summary.csv",
                  "breast_SCAN_B/data_for_summary.csv")


# GET DATA AND CREATE MASTER TABLE ---------------------------------------------------------------------------------------------------

setwd(data_directory_path)

all_datasets <- NULL

for (i in 1:length(datasets_to_use)) {
  all_datasets <- bind_rows(all_datasets, read.csv(datasets_to_use[i]))
}

# GENERAL THINGS  ----------------------------------------------------------------------------------------------------------------

full_names_plus <- c('THCA' = 'Thyroid carcinoma (THCA)', 'KIRP' = 'Kidney renal papillary cell carcinoma (KIRP)', 'PRAD' = 'Prostate adenocarcinoma (PRAD)', 
                     'LUAD' = 'Lung adenocarcinoma (LUAD)', 'KICH' = 'Kidney chromophobe (KICH)', 'KIRC' = 'Kidney renal clear cell carcinoma (KIRC)', 
                     'PCPG' = 'Pheochromocytoma and paraganglioma (PCPG)', 'LGG' = 'Brain lower grade glioma (LGG)', 'ACC' = 'Adrenocortical carcinoma (ACC)', 
                     'LIHC' = 'Liver hepatocellular carcinoma (LIHC)', 'PAAD' = 'Pancreatic adenocarcinoma (PAAD)', 'BRCA' = 'Breast invasive carcinoma (BRCA)', 
                     'MESO' = 'Mesothelioma (MESO)', 'CHOL' = 'Cholangiocarcinoma (CHOL)', 'THYM' = 'Thymoma (THYM)', 'STAD' = 'Stomach adenocarcinoma (STAD)', 
                     'LUSC' = 'Lung squamous cell carcinoma (LUSC)', 'UVM' = 'Uveal melanoma (UVM)', 'SARC' = 'Sarcoma (SARC)', 'BLCA' = 'Bladder urothelial carcinoma (BLCA)', 
                     'OV' = 'Ovarian serous cystadenocarcinoma (OV)', 'GBM' = 'Glioblastoma multiforme (GBM)', 'HNSC' = 'Head and neck squamous cell carcinoma (HNSC)', 
                     'ESCA' = 'Esophageal carcinoma (ESCA)', 'COAD' = 'Colon adenocarcinoma (COAD)', 'READ' = 'Rectum adenocarcinoma (READ)',
                     'UCEC' = 'Uterine corpus endometrial carcinoma (UCEC)', 'UCS' = 'Uterine carcinosarcoma (UCS)',
                     'CESC' = 'Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)', 
                     'TGCT' = 'Testicular germ cell tumors (TGCT)', 'SKCM' = 'Skin cutaneous melanoma (SKCM)', 'DLBC' = 'Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)',
                     'BRCA-G' = '[GOBO] Breast invasive carcinoma (BRCA)', 'BRCA-S' = '[SCAN B] Breast invasive carcinoma (BRCA)')

location_system_flip <- c("SARC", "SKCM", "DLBC", "THCA", "PCPG", "ACC", "BLCA", "KIRP", 
                          "KIRC", "KICH", "PRAD", "TGCT", "CESC", "UCS", "UCEC",
                          "OV", "BRCA - SCAN B", "BRCA - GOBO", "BRCA", "MESO", "THYM", "LUSC", "LUAD", "READ",
                          "COAD", "PAAD", "CHOL", "LIHC", "STAD", "ESCA", "HNSC", "LGG", "GBM", "UVM")

tru_percent_flip <- c("SKCM", "DLBC", "TGCT", "CESC", "UCS", "UCEC", "READ", "COAD", "ESCA", 
                   "HNSC", "GBM", "OV", "BLCA", "SARC", "UVM", "LUSC", "STAD", "THYM", 
                   "CHOL", "MESO", "BRCA", "PAAD", "LIHC", "BRCA-G", "ACC", "BRCA-S", "LGG", "PCPG",
                   "KIRC", "KICH", "LUAD", "PRAD", "KIRP", "THCA")

# BARPLOT - SAMPLE SIZES ----------------------------------------------------------------------------------------------------------------

# order by position in the body
all_datasets %>%
  mutate(cancer_order = factor(cancer_type, levels=c(location_system_flip))) %>%
  ggplot(aes(x=cancer_order, y=clams_number_samples)) +
  geom_bar(stat="identity", fill="darkgrey") +
  coord_flip() +
  labs(x ="", y = "Number of samples") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position="none") +
  scale_x_discrete(labels = full_names_plus)

# order by sample size
all_datasets %>%
  ggplot(aes(x=reorder(cancer_type, -total_samples), y=clams_number_samples)) +
  geom_bar(stat="identity", fill="darkgrey") +
  coord_flip() +
  labs(x ="", y = "Number of samples") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position="none") +
  scale_x_discrete(labels = full_names_plus)

# order by sample size - lollipop
all_datasets %>%
  ggplot(aes(x=reorder(cancer_type, -total_samples), y=total_samples)) +
  geom_segment( aes(x=reorder(cancer_type, -total_samples), xend=reorder(cancer_type, -total_samples), y=0, yend=total_samples), color="darkgrey") +
  geom_point(size=2, color="darkgrey") +
  coord_flip() +
  labs(x ="", y = "Number of samples") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position="none") +
  scale_x_discrete(labels = full_names_plus)


# BARPLOT - TRU X NONTRU PERCENTAGES ----------------------------------------------------------------------------------------------------

all_datasets %>%
  mutate(cancer_order = factor(cancer_type, levels=c(tru_percent_flip))) %>%
  ggplot(aes(x=cancer_order, y=sample_percent)) +
  geom_bar(stat="identity", aes(fill=clams_class)) +
  theme_minimal() +
  coord_flip() +
  labs(x ="", y = "Sample percentage") +
  geom_text(aes(x=cancer_type, y=105, label=total_samples), size=3, color="grey50") +
  geom_text(aes(x=34, y=113, label="= n"), size=3, color="grey50") +
  theme(axis.text.x = element_text(size = 11),
        panel.grid.major.x = element_line(),
        panel.grid.major.y = element_blank()) +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), guide=guide_legend(reverse = TRUE), name="CLAMS") +
  scale_x_discrete(labels = full_names_plus)

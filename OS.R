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
library(survival)
library(survminer)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX




### BEGINNING OF OS FROM PROLIF ---------------------

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


### END OF OS FROM PROLIF ---------------------











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

#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO RUN COX REGRESSION ANALYSIS
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
# (1B) OS FROM ANNOTATION
# CLAMS RESULTS
# PROLIFERATION RESULTS
# ROR RESULTS
# PAM50 RESULTS


patient_os_all_datasets <- patient_os_all_datasets %>% mutate(groups.to.analyze = paste(cancer.type, dataset, sep="_"))
patient_os_all_datasets$groups.to.analyze <- factor(patient_os_all_datasets$groups.to.analyze)
groups.to.analyze <- levels(patient_os_all_datasets$groups.to.analyze)


# GET ALL TIME INFO TO SAME TYPE  (M/Y)  --------------------------

# mutate patient OS info to months and years
patient_os_all_datasets$OS.time.type <- as.character(patient_os_all_datasets$OS.time.type)
patient_os_all_datasets <- patient_os_all_datasets %>%
  mutate(OS.time.years = case_when(OS.time.type == 'days' ~ (OS.time / 365.2425),
                                   OS.time.type == 'months' ~ (OS.time / 30.436875),
                                   OS.time.type == 'years' ~ OS.time,
                                   TRUE ~ OS.time),
         OS.time.months = case_when(OS.time.type == 'days' ~ (OS.time / 30.436875),
                                    OS.time.type == 'months' ~ OS.time,
                                    OS.time.type == 'years' ~ (OS.time * 12),
                                    TRUE ~ OS.time))


# CENSOR AT SPECIFIC TIME POINT  ----------------------------------

censor_data <- function (time_to_censor, unit_to_censor) {
  # get new names and what columns to use for censoring
  time_column_to_use <- case_when(unit_to_censor == "y" ~ "OS.time.years",
                                  unit_to_censor == "m" ~ "OS.time.months")
  censored_time_column <- case_when(unit_to_censor == "y" ~ paste0("OS.time.years.", time_to_censor, unit_to_censor),
                                    unit_to_censor == "m" ~ paste0("OS.time.months.", time_to_censor, unit_to_censor))
  censored_event_column <- paste0("OS.event.", time_to_censor, unit_to_censor)
  
  # censor
  patient_os_all_datasets <- patient_os_all_datasets %>%
    # if OS number > time_to_censor, set time to cap, else keep original number
    mutate(censored_time =
             ifelse(get(time_column_to_use) > as.numeric(time_to_censor),
                    as.numeric(time_to_censor), get(time_column_to_use)),
           # if OS number > time_to_censor, set event to 0, else keep original event
           censored_event = 
             ifelse(get(time_column_to_use) > as.numeric(time_to_censor),
                    0, OS.event))
  names(patient_os_all_datasets)[names(patient_os_all_datasets) == "censored_time"] <- censored_time_column
  names(patient_os_all_datasets)[names(patient_os_all_datasets) == "censored_event"] <- censored_event_column  
  return(patient_os_all_datasets)
}               

patient_os_all_datasets <- censor_data(5, "y")
patient_os_all_datasets <- censor_data(10, "y")


# CLAMS  -------

setwd("/media/deboraholi/Data/LUND/9 THESIS/0_clams/Cox/")

patients_os_clams <- left_join(patient_os_all_datasets, patient_annotation_clams,
                               by = c("sample.id", "cancer.type", "dataset"))

groups.to.analyze.clams <- subset(patients_os_clams, clams.class == "TRU") %>% pull(groups.to.analyze) %>% unique() %>% as.character()

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

# sink("clams.txt") - to save output to a file
# individual Cox forest plots
for (type_to_analyze in groups.to.analyze.clams) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams, groups.to.analyze == type_to_analyze)
  type_subset$clams.class <- factor(type_subset$clams.class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~clams.class, data = type_subset)))
  forest <- ggforest(coxph(formula = overall_surv_object~clams.class, data = type_subset), 
                 main = paste("Hazard ratio", type_to_analyze))
  print(forest)
  current_filename <- paste0("CLAMS Cox ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(forest), width=6, height=3, dpi=300)
}
# sink()

# PROLIFERATION  -------

# 25th percentile
setwd("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/25th percentile/Cox/")

patients_os_prolif <- left_join(patient_os_all_datasets, patient_annotation_prolif,
                                by = c("sample.id", "cancer.type", "dataset", "groups.to.analyze"))
patients_os_prolif$groups.to.analyze <- factor(patients_os_prolif$groups.to.analyze)
groups.to.analyze <- levels(factor(patients_os_prolif$groups.to.analyze))

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

# sink("prolif.txt") - to save output to a file
# individual Cox forest plots
for (type_to_analyze in groups.to.analyze) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_prolif, groups.to.analyze == type_to_analyze)
  type_subset$prolif.group.karl <- factor(type_subset$prolif.group.karl, levels=c("Low", "High"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~prolif.group.karl, data = type_subset)))
  forest <- ggforest(coxph(formula = overall_surv_object~prolif.group.karl, data = type_subset), 
                     main = paste("Hazard ratio", type_to_analyze))
  print(forest)
  current_filename <- paste0("Prolif Cox ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(forest), width=6, height=3, dpi=300)
}
# sink()


# ROR  -------

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/Cox/")

patients_os_ror <- left_join(patient_os_all_datasets, patient_annotation_ror,
                             by = c("sample.id", "cancer.type", "dataset"))

# extract only the ones with more than one ROR class
groups.to.analyze.ror <- NULL
for (type_to_analyze in levels(factor(patients_os_ror$groups.to.analyze))) {
  number_of_classes <- patients_os_ror %>% filter(groups.to.analyze == type_to_analyze) %>% 
                            group_by(ror.red.hl.class) %>% tally() %>% nrow()
  if (number_of_classes > 1) {
    groups.to.analyze.ror <- c(groups.to.analyze.ror, type_to_analyze)
  }
}

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"
groups.to.analyze.clams <- subset(patients_os_clams, clams.class == "TRU") %>% pull(groups.to.analyze) %>% unique() %>% as.character()
# sink("ror.txt") - to save output to a file
# individual Cox forest plots
for (type_to_analyze in groups.to.analyze.ror) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_ror, groups.to.analyze == type_to_analyze)
  classes <- type_subset %>% pull(ror.red.hl.class)%>% unique() %>% as.character()
  type_subset$ror.red.hl.class <- factor(type_subset$ror.red.hl.class, levels=c("Low", "Medium", "High"))
  # if Low not in sample, relevel with Medium
  if (!"Low" %in% classes) {
    type_subset$ror.red.hl.class <- factor(type_subset$ror.red.hl.class, levels=c("Medium", "High"))
  }
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~ror.red.hl.class, data = type_subset)))
  forest <- ggforest(coxph(formula = overall_surv_object~ror.red.hl.class, data = type_subset), 
                     main = paste("Hazard ratio", type_to_analyze))
  print(forest)
  current_filename <- paste0("ROR Cox ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(forest), width=6, height=3, dpi=300)
}
# sink()


# PAM50  -------

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/Cox/")

patients_os_pam50 <- left_join(patient_os_all_datasets, patient_annotation_pam50,
                                by = c("sample.id", "cancer.type", "dataset"))

# extract only the ones with more than one PAM50 class - all
groups.to.analyze.pam50 <- NULL
for (type_to_analyze in levels(factor(patients_os_pam50$groups.to.analyze))) {
  number_of_classes <- patients_os_pam50 %>% filter(groups.to.analyze == type_to_analyze) %>% 
    group_by(pam50.red.class) %>% tally() %>% nrow()
  if (number_of_classes > 1) {
    groups.to.analyze.pam50 <- c(groups.to.analyze.pam50, type_to_analyze)
  }
}

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

# sink("pam50.txt") - to save output to a file
# individual Cox forest plots
for (type_to_analyze in groups.to.analyze.pam50) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_pam50, groups.to.analyze == type_to_analyze)
  type_subset$pam50.red.class <- factor(type_subset$pam50.red.class, levels=c("LumA", "LumB", "Her2", "Basal"))
  classes <- type_subset %>% pull(pam50.red.class)%>% unique() %>% as.character()
  # relevel if needed
  if (!"LumA" %in% classes & !"LumB" %in% classes) {
    type_subset$pam50.red.class <- factor(type_subset$pam50.red.class, levels=c("Her2", "Basal"))
  }
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~pam50.red.class, data = type_subset)))
  forest <- ggforest(coxph(formula = overall_surv_object~pam50.red.class, data = type_subset), 
                     main = paste("Hazard ratio", type_to_analyze))
  print(forest)
  current_filename <- paste0("PAM50 Cox ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(forest), width=6, height=3, dpi=300)
}
# sink()
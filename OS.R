#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# SCRIPT TO CENSOR OS DATA AND PLOT
#
#

# CLEAR ENVIRONMENT ------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ----------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(viridisLite)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1B) OS FROM ANNOTATION - always
# if needed:
# CLAMS RESULTS
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


# OVERALL SURVIVAL CLAMS  -------------------------------------------------

setwd("/media/deboraholi/Data/LUND/9 THESIS/0_clams/OS 5y/")

patients_os_clams <- left_join(patient_os_all_datasets, patient_annotation_clams,
                               by = c("sample.id", "cancer.type", "dataset"))

groups.to.analyze.clams <- subset(patients_os_clams, clams.class == "TRU") %>% pull(groups.to.analyze) %>% unique() %>% as.character()

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

for (type_to_analyze in groups.to.analyze.clams) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams, groups.to.analyze == type_to_analyze)
  type_subset$clams.class <- factor(type_subset$clams.class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  # separating estimates by TRU or NonTRU
  fit_clams <- survfit(overall_surv_object~clams.class, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column]], 
                      type_subset[[event_column]])~type_subset$clams.class,))
  
  # plot
  # simple for comparing info and making sure names are correct
  # print(ggsurvplot(fit_clams, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  current_plot <- ggsurvplot(fit_clams, data=type_subset,
                   palette=c("darkorange1", "deepskyblue3"),
                   title=paste0("Overall Survival (", type_to_analyze, ")"), 
                   xlab=paste0("Time (", time_type, ")"),
                   censor.shape=124, censor.size=3,
                   pval=TRUE, pval.coord=c(0,0.1),
                   surv.median.line="hv",
                   risk.table=TRUE,
                   risk.table.fontsize = 4,
                   tables.theme = theme_survminer(font.main = 14),
                   legend="none", legend.title="CLAMS",
                   legend.labs=c("TRU", "NonTRU"))
  print(current_plot)
  current_filename <- paste0("OS 5y ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}


setwd("/media/deboraholi/Data/LUND/9 THESIS/0_clams/OS 10y BRCA/")

groups.brca <- subset(patients_os_clams, cancer.type == "BRCA") %>% pull(groups.to.analyze) %>% unique() %>% as.character()

time_column <- "OS.time.years.10y"
event_column <- "OS.event.10y"
time_type <- "years"

for (type_to_analyze in groups.brca) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams, groups.to.analyze == type_to_analyze)
  type_subset$clams.class <- factor(type_subset$clams.class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  # separating estimates by TRU or NonTRU
  fit_clams <- survfit(overall_surv_object~clams.class, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column]], 
                      type_subset[[event_column]])~type_subset$clams.class,))
  
  # plot
  # simple for comparing info and making sure names are correct
  # print(ggsurvplot(fit_clams, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  current_plot <- ggsurvplot(fit_clams, data=type_subset,
                             palette=c("darkorange1", "deepskyblue3"),
                             title=paste0("Overall Survival (", type_to_analyze, ")"), 
                             xlab=paste0("Time (", time_type, ")"),
                             censor.shape=124, censor.size=3,
                             pval=TRUE, pval.coord=c(0,0.1),
                             surv.median.line="hv",
                             risk.table=TRUE,
                             risk.table.fontsize = 4,
                             tables.theme = theme_survminer(font.main = 14),
                             legend="none", legend.title="CLAMS",
                             legend.labs=c("TRU", "NonTRU"))
  print(current_plot)
  current_filename <- paste0("OS 10y ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}


# PROLIFERATION  -------------------------------------------------

setwd("/media/deboraholi/Data/LUND/9 THESIS/1_proliferation/Karlsson/25th percentile/OS 5y/")

patients_os_prolif <- left_join(patient_os_all_datasets, patient_annotation_prolif,
                               by = c("sample.id", "cancer.type", "dataset", "groups.to.analyze"))

groups.to.analyze <- levels(factor(patients_os_prolif$groups.to.analyze))

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"


for (a.group in groups.to.analyze) {
  print(a.group)
  group_subset <- subset(patients_os_prolif, groups.to.analyze == a.group)
  group_subset$prolif.group.karl <- factor(group_subset$prolif.group.karl, levels=c("Low", "High"))
  # make the object
  overall_surv_object <- Surv(time=group_subset[[time_column]], 
                              event=group_subset[[event_column]])
  # separating estimates by Low or High
  fit_prolif <- survfit(overall_surv_object~prolif.group.karl, data=group_subset)
  # testing by proliferation group
  print(survdiff(Surv(group_subset[[time_column]], 
                      group_subset[[event_column]])~group_subset$prolif.group.karl,))
  # plot
  # simple for comparing info and making sure names are correct
  # print(ggsurvplot(fit_prolif, data=group_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  current_plot <- ggsurvplot(fit_prolif, data=group_subset,
                             palette=c("#FDE725", "#26828E"),
                             title=paste0("Overall Survival (", a.group, ")"), 
                             xlab=paste0("Time (", time_type, ")"),
                             censor.shape=124, censor.size=3,
                             pval=TRUE, pval.coord=c(0,0.1),
                             surv.median.line="hv",
                             risk.table=TRUE,
                             risk.table.fontsize = 4,
                             tables.theme = theme_survminer(font.main = 14),
                             legend="none", legend.title="Proliferation",
                             legend.labs=c("Low", "High"))
  print(current_plot)
  current_filename <- paste0("OS 5y ", a.group, ".png")
  ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}


# ROR  -------------------------------------------------
# but do with relapse free instead of OS?

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/ROR/OS 5y/")

patients_os_ror <- left_join(patient_os_all_datasets, patient_annotation_ror,
                                by = c("sample.id", "cancer.type", "dataset"))

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

# ROR reduced set
plot_ror_red_os <- function(ror_group, ror_labels, ror_colors) {
  for (a.group in ror_group) {
    print(a.group)
    group_subset <- subset(patients_os_ror, groups.to.analyze == a.group)
    group_subset$ror.red.hl.class <- factor(group_subset$ror.red.hl.class, levels=ror_labels)
    # make the object
    overall_surv_object <<- Surv(time=group_subset[[time_column]], # bad practice, but fixed my problem
                                event=group_subset[[event_column]])
    # separating estimates by Low, Medium or High
    fit_prolif <- survfit(overall_surv_object~ror.red.hl.class, data=group_subset)
    # testing by ror group
    print(survdiff(Surv(group_subset[[time_column]], 
                          group_subset[[event_column]])~group_subset$ror.red.hl.class,))
    # plot
    # simple for comparing info and making sure names are correct
    # print(ggsurvplot(fit_prolif, data=group_subset, title="Overall Survival", risk.table=TRUE))
    # good looking one
    current_plot <- ggsurvplot(fit_prolif, data=group_subset,
                               palette=ror_colors,
                               title=paste0("Overall Survival (", a.group, ")"),
                               xlab=paste0("Time (", time_type, ")"),
                               censor.shape=124, censor.size=3,
                               pval=TRUE, pval.coord=c(0,0.1),
                               surv.median.line="hv",
                               risk.table=TRUE,
                               risk.table.fontsize = 4,
                               tables.theme = theme_survminer(font.main = 14),
                               legend="none", legend.title="ROR",
                               legend.labs=ror_labels)
    print(current_plot)
    current_filename <- paste0("OS 5y ", a.group, ".png")
    ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
  }
}

ror.groups.lmh <- NULL
ror.groups.lm <- NULL
ror.groups.mh <- NULL
ror.groups.hl <- NULL
ror.groups.one <- NULL

# divide plots by ROR classes in data
for (a.group in groups.to.analyze) {
  ror_classes <- subset(patients_os_ror, groups.to.analyze == a.group) %>% pull(ror.red.hl.class) %>% unique() %>% as.character()
  if (length(ror_classes) == 3) {
    ror.groups.lmh <- c(ror.groups.lmh, a.group)
  }
  else if (length(ror_classes) == 2) {
    if ("Low" %in% ror_classes & "Medium" %in% ror_classes) {
      ror.groups.lm <- c(ror.groups.lm, a.group)
    } else if ("Medium" %in% ror_classes & "High" %in% ror_classes) {
      ror.groups.mh <- c(ror.groups.mh, a.group)
    } else if ("Low" %in% ror_classes & "High" %in% ror_classes) {
      ror.groups.hl <- c(ror.groups.hl, a.group)
    }
  }
  else if (length(ror_classes) == 1) {
    ror.groups.one <- c(ror.groups.one, a.group)
  }
}

ror_labels <- c("Low", "Medium", "High")
ror_colors <- c("#bcb9b9", "#6a6868", "#030303")

plot_ror_red_os(ror.groups.lmh, ror_labels, ror_colors)
plot_ror_red_os(ror.groups.lm, ror_labels[1:2], ror_colors[1:2])
plot_ror_red_os(ror.groups.mh, ror_labels[2:3], ror_colors[2:3])
# plot_ror_red_os(ror.groups.hl, "Low", "#bcb9b9") # empty


# PAM50  --------------------------------------------------

setwd("/media/deboraholi/Data/LUND/9 THESIS/3_brca_ssps/PAM50/OS 5y/")

patients_os_pam50 <- left_join(patient_os_all_datasets, patient_annotation_pam50,
                             by = c("sample.id", "cancer.type", "dataset"))
patients_os_pam50$pam50.red.class <- factor(patients_os_pam50$pam50.red.class, 
                                            levels=c("LumA", "LumB", "Her2", "Basal"))

# BRCA only
brca_subset_os_pam50 <- subset(patients_os_pam50, cancer.type == "BRCA")

groups.to.analyze <- levels(factor(brca_subset_os_pam50$groups.to.analyze))

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

# reduced
for (a.group in groups.to.analyze) {
  print(a.group)
  group_subset <- subset(brca_subset_os_pam50, groups.to.analyze == a.group)
  # make the object
  overall_surv_object <- Surv(time=group_subset[[time_column]], 
                              event=group_subset[[event_column]])
  # separating estimates by PAM50
  fit_prolif <- survfit(overall_surv_object~pam50.red.class, data=group_subset)
  # testing by proliferation group
  print(survdiff(Surv(group_subset[[time_column]], 
                      group_subset[[event_column]])~group_subset$pam50.red.class,))
  # plot
  # simple for comparing info and making sure names are correct
  # print(ggsurvplot(fit_prolif, data=group_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  current_plot <- ggsurvplot(fit_prolif, data=group_subset,
                             palette=plasma(4, begin = 0.3, direction = -1),
                             title=paste0("Overall Survival (", a.group, ")"),
                             xlab=paste0("Time (", time_type, ")"),
                             censor.shape=124, censor.size=3,
                             pval=TRUE, pval.coord=c(0,0.1),
                             surv.median.line="hv",
                             risk.table=TRUE,
                             risk.table.fontsize = 4,
                             tables.theme = theme_survminer(font.main = 14),
                             legend="none", legend.title="PAM50",
                             legend.labs=c("LumA", "LumB", "HER2", "Basal"))
  print(current_plot)
  current_filename <- paste0("OS 5y ", a.group, ".png")
  ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}
























### BEGINNING OF COX FROM PROLIF ---------------------


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

### END OF COX FROM PROLIF ---------------------




















# Hazard ratios


setwd("/media/deboraholi/Data/LUND/9 THESIS/data")
patient_annotation_scanb <- readRDS("breast_SCAN_B/patient_clams.rds")
patient_annotation_gobo <- readRDS("breast_GOBO/patient_clams.rds")
patient_annotation_tcga <- readRDS("all_tcga/patient_clams.rds")
#patient_annotation_tcga <- readRDS("all_tcga/patient_clams.rds")
patient_annotation_tcga_brca <- subset(patient_annotation_tcga, cancer.type == "BRCA")
patient_annotation_tcga <- patient_annotation_tcga_brca

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

# FOR ONLY ONE CANCER TYPE IN THE DATASET
# analysis
patient_annotation$clams_class <- factor(patient_annotation$clams_class, levels=c("TRU", "NonTRU"))
# make the object
overall_surv_object <- Surv(time=patient_annotation[[time_column_to_use]], 
                            event=patient_annotation[[event_column_to_use]])
# separating estimates by TRU or NonTRU
fit_clams <- survfit(overall_surv_object~clams_class, data=patient_annotation)
# testing by clams class
print(survdiff(Surv(patient_annotation[[time_column_to_use]], 
                    patient_annotation[[event_column_to_use]])~patient_annotation$clams_class,))

# plot
# simple for comparing info and making sure names are correct
print(ggsurvplot(fit_clams, data=patient_annotation, title="Overall Survival", risk.table=TRUE))
# good looking one
print(ggsurvplot(fit_clams, data=patient_annotation,
                 palette=c("darkorange1", "deepskyblue3"),
                 title=paste0("Overall Survival (", config$overall_survival_name, ")"), 
                 xlab=paste0("Time (", config$patient_annotation_info$overall_survival_time_column_unit_desired, ")"),
                 censor.shape=124, censor.size=3,
                 pval=TRUE, pval.coord=c(0,0.1),
                 surv.median.line="hv",
                 risk.table=TRUE,
                 risk.table.fontsize = 4,
                 tables.theme = theme_survminer(font.main = 14),
                 legend="none", legend.title="CLAMS",
                 legend.labs=c("TRU", "NonTRU")))
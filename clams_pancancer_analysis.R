#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO PREPARE CANCER DATA, RUN CLAMS PREDICTOR, ANALYSE DATA
#
#


# CLEAR ENVIRONMENT ----------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(yaml)



# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# don't source or it crashes...
# (1) PATIENT ANNOTATION # only for TCGA
# (1A) ONLY PATIENT IDS FROM ANNOTATION
# (2) GENE TABLE INFORMATION
# (3) GEX



# CLAMS ANALYSIS  -----------------------------------

library(CLAMS)
library(e1071)


# GOBO

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_gobo.yml")

# Create a gene symbol list from the reporter ids present in the gex_matrix  
gene_symbols_gobo <- gex_matrix_gobo[,1:2]
gene_symbols_gobo <- rownames_to_column(data.frame(gene_symbols_gobo), var = "reporterId")
gene_symbols_gobo <- left_join(gene_symbols_gobo, gene_table_gobo[,c("reporterId","geneSymbol")])
gene_symbols_gobo <- gene_symbols_gobo$geneSymbol

rm(gene_table_gobo)

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_gobo) <- gene_symbols_gobo

# run CLAMS predictor
clams_result_gobo <- applyCLAMS(gex_matrix_gobo, gene_symbols_gobo)
# extract the classification result from CLAMS
clams.class <- clams_result_gobo$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_gobo <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_gobo)
rm(gene_symbols_gobo)



# SCAN-B

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_brca_scanb.yml")

# Create a gene symbol list from gex_matrix
gene_symbols_scanb <- gex_matrix_scanb[,1:2]
gene_symbols_scanb <- rownames_to_column(data.frame(gex_matrix_scanb), var = "geneSymbol")
gene_symbols_scanb <- gene_symbols_scanb$geneSymbol

# run CLAMS predictor
clams_result_scanb <- applyCLAMS(gex_matrix_scanb, gene_symbols_scanb)
# extract the classification result from CLAMS
clams.class <- clams_result_scanb$cl[,]
clams.class <- data.frame(as.list(clams.class))
clams.class <- t(clams.class)
clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")

clams_class_scanb <- inner_join(patient_ids_all_datasets, clams.class, by="sample.id")

rm(gex_matrix_scanb)
rm(gene_symbols_scanb)


# TCGA

setwd("/media/deboraholi/Data/LUND/9 THESIS/src")
config <- yaml.load_file("config_alltcga.yml")

# Create a gene symbol list from the reporter ids present in the gex_matrix  
gene_symbols_tcga <- gex_matrix_tcga[,1:2]
gene_symbols_tcga <- rownames_to_column(data.frame(gene_symbols_tcga), var = "reporterId")
gene_symbols_tcga <- left_join(gene_symbols_tcga, gene_table_tcga[,c("ENSG", "SYMBOL")], by=c("reporterId" = "ENSG"))
gene_symbols_tcga <- gene_symbols_tcga$SYMBOL

rm(gene_table_tcga)

# Replace VEGFD for FIGF
gene_symbols_tcga <- as.list(gene_symbols_tcga)
gene_symbols_tcga <- rapply(gene_symbols_tcga, function(x) ifelse(x=="VEGFD","FIGF",x), how="replace")

# Rename gex_matrix rows with gene symbols instead of reporter ids
row.names(gex_matrix_tcga) <- gene_symbols_tcga

rm(gene_symbols_tcga)

# run CLAMS predictor for subsets of each cancer type
cancer_types <- levels(factor(patient_annotation_tcga$cancer.type))
for (cancer_type in cancer_types) {
  
  # Subset data according to cancer type
  print(cancer_type)
  sample_list <- subset(patient_annotation_tcga, cancer.type == cancer_type) %>% pull(sample_id)
  col_num <- which(colnames(gex_matrix_tcga) %in% sample_list)
  new_gex_matrix <- gex_matrix_tcga[,c(col_num)]
  gene_symbols <- row.names(new_gex_matrix)
  
  # run CLAMS predictor
  clams_result <- applyCLAMS(new_gex_matrix, gene_symbols)
  # extract the classification result from CLAMS
  clams.class <- clams_result$cl[,]
  
  clams.class <- data.frame(as.list(clams.class))
  clams.class <- t(clams.class)
  clams.class <- rownames_to_column(data.frame(clams.class), var = "sample.id")
  
  if (exists("clams_class_tcga")) {
    clams_class_tcga <- bind_rows(clams_class_tcga, clams.class) } else {
      clams_class_tcga <- clams.class
    }
  
  rm(new_gex_matrix)
}

clams_class_tcga <- inner_join(patient_ids_all_datasets, clams_class_tcga, by="sample.id")

rm(patient_annotation_tcga)
rm(gex_matrix_tcga)


# combine all datasets
patient_annotation_clams <- bind_rows(clams_class_gobo, clams_class_scanb, clams_class_tcga)
# transmute bronchioid to TRU and nonbronchioid to NonTRU
patient_annotation_clams <- patient_annotation_clams %>% 
                              mutate(clams.class = case_when(clams.class == "bronchioid" ~ "TRU",
                                     clams.class == "nonbronchioid" ~ "NonTRU",
                                     TRUE ~ clams.class))
patient_annotation_clams$clams.class <- factor(patient_annotation_clams$clams.class, levels=c("TRU", "NonTRU"))


# RESULT
# Save final table
write.csv(patient_annotation_clams, "/media/deboraholi/Data/LUND/9 THESIS/0_clams/clams_all_samples.csv", row.names=FALSE)


  
# CHANGE FROM HERE ------------------



# GENERAL VISUALIZATION OF PATIENT ANNOTATION INFORMATION BY CLAMS CLASS (ONE TYPE IN DATASET) --------------------------------------------------

# Barplots
if (length(config$see_variables_by_clams$barplot) > 0) {
  for (i in 1:length(config$see_variables_by_clams$barplot)) {
    variable <- config$see_variables_by_clams$barplot[i]
    this_plot <- ggplot(data=subset(patient_annotation, !is.na(get(variable))), aes_string(x=variable)) +
      geom_bar(aes(fill=clams_class)) +
      labs(y ="Number of samples") +
      theme_minimal() +
      scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}


# Histograms
if (length(config$see_variables_by_clams$histogram) > 0) {
  for (i in 1:length(config$see_variables_by_clams$histogram)) {
    variable <- config$see_variables_by_clams$histogram[i]
    this_plot <- ggplot(data=patient_annotation, aes_string(x=variable)) +
      geom_histogram(aes(fill=clams_class)) +
      labs(y ="Number of samples") +
      theme_minimal() +
      scale_fill_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}


# Frequency polygons
if (length(config$see_variables_by_clams$histogram) > 0) {
  for (i in 1:length(config$see_variables_by_clams$histogram)) {
    variable <- config$see_variables_by_clams$histogram[i]
    this_plot <- ggplot(data=patient_annotation, aes(get(variable), stat(density), colour=clams_class)) +
      geom_freqpoly(bins=15) +
      labs(x = variable) +
      theme_minimal() +
      scale_colour_manual(values=c("deepskyblue3", "darkorange1"), name="CLAMS")
    print(this_plot)
  }
}


# SURVIVAL ANALYSIS (GENERAL) ----------------------------------------------------------------------------------------------------------------------

# DATA PREP

# MUTATE PATIENT_ANNOTATION SO OS TIME IS IN MONTHS OR YEARS INSTEAD OF DAYS IF NEEDED
if (config$patient_annotation_info$overall_survival_time_column_unit_table == config$patient_annotation_info$overall_survival_time_column_unit_desired) {
  os_time_column_to_use <- config$patient_annotation_info$overall_survival_time_column_name
} else if (config$patient_annotation_info$overall_survival_time_column_unit_table == 'days' & 
    config$patient_annotation_info$overall_survival_time_column_unit_desired == 'years') {
  patient_annotation <- mutate(patient_annotation, os_time_years = get(config$patient_annotation_info$overall_survival_time_column_name) / 365.2425)
  new_os_time_column_name <- paste0(config$patient_annotation_info$overall_survival_time_column_name, "_years")
  names(patient_annotation)[names(patient_annotation) == "os_time_years"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else if (config$patient_annotation_info$overall_survival_time_column_unit_table == 'months' & 
    config$patient_annotation_info$overall_survival_time_column_unit_desired == 'years') {
  patient_annotation <- mutate(patient_annotation, os_time_years = get(config$patient_annotation_info$overall_survival_time_column_name) / 12)
  new_os_time_column_name <- paste0(config$patient_annotation_info$overall_survival_time_column_name, "_years")
  names(patient_annotation)[names(patient_annotation) == "os_time_years"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else if (config$patient_annotation_info$overall_survival_time_column_unit_table == 'days' & 
         config$patient_annotation_info$overall_survival_time_column_unit_desired == 'months') {
  patient_annotation <- mutate(patient_annotation, os_time_months = get(config$patient_annotation_info$overall_survival_time_column_name) / 30.436875)
  new_os_time_column_name <- paste0(config$patient_annotation_info$overall_survival_time_column_name, "_months")
  names(patient_annotation)[names(patient_annotation) == "os_time_months"] <- new_os_time_column_name
  os_time_column_to_use <- new_os_time_column_name
} else {
  print('Something went wrong when changing the time unit')
}


# CENSOR DATA AT SPECIFIC NUMBER OF DAYS/MONTHS/YEARS
if (!config$survival_analysis$censor_at_timepoint == 'no' & config$survival_analysis$censor_at_timepoint > 0) {
  number_to_cap_at <- config$survival_analysis$censor_at_timepoint
  # if OS number > number_to_cap_at, set time to cap, else keep original number
  patient_annotation <- mutate(patient_annotation, new_os_time_column =
                                 ifelse(
                                   get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                   as.numeric(number_to_cap_at), get(os_time_column_to_use)))
  # if OS number > number_to_cap_at, set event to 0, else keep original event
  patient_annotation <- mutate(patient_annotation, new_os_event_column =
                                 ifelse(
                                   get(os_time_column_to_use) > as.numeric(number_to_cap_at),
                                   0, get(config$patient_annotation_info$overall_survival_event_column_name)))
  new_os_time_column_name <- paste(os_time_column_to_use,
                                   config$survival_analysis$censor_at_timepoint, sep="_")
  new_os_event_column_name <- paste(config$patient_annotation_info$overall_survival_event_column_name,
                                    config$survival_analysis$censor_at_timepoint, sep="_")
  names(patient_annotation)[names(patient_annotation) == "new_os_time_column"] <- new_os_time_column_name
  names(patient_annotation)[names(patient_annotation) == "new_os_event_column"] <- new_os_event_column_name  
}



# ANALYSIS

library(survival)
library(survminer)

# S1
# Overall Survival curve with the option of time censoring

# GET THE CORRECT COLUMN FOR OS TIME (DAYS/MONTHS/YEARS, CENSORED/NOT)
if (!config$survival_analysis$censor_at_timepoint == 'no' & config$survival_analysis$censor_at_timepoint > 0) {
  time_column_to_use <- new_os_time_column_name
  event_column_to_use <- new_os_event_column_name
} else if (config$survival_analysis$censor_at_timepoint == 'no') {
  time_column_to_use <- os_time_column_to_use
  event_column_to_use <- config$patient_annotation_info$overall_survival_event_column_name
} else {
  print("You did not specify a valid time to use for survival data. Please review.")
}

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


# FOR MORE THAN ONE CANCER TYPE IN THE SAME DATASET
# analysis
for (type_to_analyze in config$survival_analysis$cancer_types_TRU_OS) {
  type_subset <- subset(patient_annotation, cancer.type == type_to_analyze)
  type_subset$clams_class <- factor(type_subset$clams_class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column_to_use]], 
                              event=type_subset[[event_column_to_use]])
  # separating estimates by TRU or NonTRU
  fit_clams <- survfit(overall_surv_object~clams_class, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column_to_use]], 
                      type_subset[[event_column_to_use]])~type_subset$clams_class,))
  
  # plot
  # simple for comparing info and making sure names are correct
  print(ggsurvplot(fit_clams, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  print(ggsurvplot(fit_clams, data=type_subset,
                   palette=c("darkorange1", "deepskyblue3"),
                   title=paste0("Overall Survival (", type_to_analyze, ")"), 
                   xlab=paste0("Time (", config$patient_annotation_info$overall_survival_time_column_unit_desired, ")"),
                   censor.shape=124, censor.size=3,
                   pval=TRUE, pval.coord=c(0,0.1),
                   surv.median.line="hv",
                   risk.table=TRUE,
                   risk.table.fontsize = 4,
                   tables.theme = theme_survminer(font.main = 14),
                   legend="none", legend.title="CLAMS",
                   legend.labs=c("TRU", "NonTRU")))
}




# SURVIVAL ANALYSIS (SPECIFIC VARIABLES) ------------------------------------------------------------------------------------------------------------------
# Overall Survival curve for specific variables with the option of time censoring

# GET THE CORRECT COLUMN FOR OS TIME (CENSORED OR NOT)
# Assumes you already ran the block of converting time units AND censoring 
if (length(config$survival_after_clams) > 0) {
  if (!config$survival_analysis$censor_at_timepoint == 'no' & config$survival_analysis$censor_at_timepoint > 0) {
    time_column_to_use <- new_os_time_column_name
    event_column_to_use <- new_os_event_column_name
  } else if (config$survival_analysis$censor_at_timepoint == 'no') {
    time_column_to_use <- os_time_column_to_use
    event_column_to_use <- config$patient_annotation_info$overall_survival_event_column_name
  } else {
    print("You did not specify a valid time to use for survival data. Please review.")
  }
  patient_annotation$clams_class <- factor(patient_annotation$clams_class, levels=c("TRU", "NonTRU"))
  for (variable_category in config$survival_after_clams) {
    column_name <- variable_category[1]
    wanted_category <- variable_category[2]
    # retrieve specific data from table
    data_to_use <- subset(patient_annotation, get(column_name) == wanted_category, 
                           select=c(time_column_to_use, event_column_to_use, "clams_class"))
    # survival
    # make the object
    object_name <- paste(column_name, wanted_category, "surv_object", sep = "_")
    assign(object_name, Surv(time=data_to_use[[time_column_to_use]],
                                event=data_to_use[[event_column_to_use]]))
    # separating estimates by TRU or NonTRU
    fit_name <- paste("fit", column_name, wanted_category, sep = "_")
    assign(fit_name, survfit(get(object_name)~clams_class, data=data_to_use))
    # testing by clams class
    print(survdiff(Surv(time=data_to_use[[time_column_to_use]],
                        event=data_to_use[[event_column_to_use]])~data_to_use$clams_class,))
    
    # plot
    # simple for comparing info and making sure names are correct
    print(ggsurvplot(get(fit_name), data=data_to_use, 
                     title=paste("Overall Survival within", column_name, "-", wanted_category),
                     risk.table=TRUE))
    # good looking one
    print(ggsurvplot(get(fit_name), data=data_to_use,
                     palette=c("darkorange1", "deepskyblue3"),
                     title=paste0("Overall Survival (", column_name, " - ", wanted_category, ")"),
                     xlab=paste0("Time (", config$patient_annotation_info$overall_survival_time_column_unit_desired, ")"),
                     censor.shape=124, censor.size=3,
                     pval=TRUE, pval.coord=c(0,0.1),
                     surv.median.line="hv",
                     risk.table=TRUE,
                     risk.table.fontsize = 4,
                     tables.theme = theme_survminer(font.main = 14),
                     legend="none", legend.title="CLAMS",
                     legend.labs=c("TRU", "NonTRU")))
  }
}


# COX ANALYSIS ------------------------------------------------------------------------------------------------------------------------------------
# Cox model with the option of time censoring ------ NOT YET WITH OPTION

# GET THE CORRECT COLUMN FOR OS TIME (CENSORED OR NOT)
if (!config$survival_analysis$censor_at_timepoint == 'no' & config$survival_analysis$censor_at_timepoint > 0) {
  time_column_to_use <- new_os_time_column_name
  event_column_to_use <- new_os_event_column_name
} else if (config$survival_analysis$censor_at_timepoint == 'no') {
  time_column_to_use <- os_time_column_to_use
  event_column_to_use <- config$patient_annotation_info$overall_survival_event_column_name
} else {
  print("You did not specify a valid time to use for survival data. Please review.")
}

# Univariate for clams_class with multiple cancer types - working
for (type_to_analyze in config$cox_regression$cancer_types) {
  print(type_to_analyze)
  type_subset <- subset(patient_annotation, cancer.type == type_to_analyze)
  type_subset$clams_class <- factor(type_subset$clams_class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column_to_use]], 
                              event=type_subset[[event_column_to_use]])
  print(summary(coxph(formula = overall_surv_object~clams_class, data = type_subset)))
  print(ggforest(coxph(formula = overall_surv_object~clams_class, data = type_subset), main = paste("Hazard ratio", type_to_analyze)))
}






  # Univariate - who knows if working now
if (config$cox_regression$do_univariate == "yes" & length(config$cox_regression$univariate_variables > 0)) {
  for (i in 1:length(config$cox_regression$univariate_variables)) {
    uni_cox_formula <- as.formula(paste("Surv(time=time_column_to_use, 
                                 event=event_column_to_use)~", config$cox_regression$univariate_variables[i]))
    uni_cox_fit <- coxph(uni_cox_formula, data=patient_annotation)
    print(uni_cox_fit)
    print(ggforest(uni_cox_fit))
    rm(uni_cox_fit)
  }
} else if (config$cox_regression$do_overall_clams == "no") {
  print("You do not want any univariate Cox regression analysis, so skipping this")
}

# Multivariate - who knows if working now
if (config$cox_regression$do_multivariate == "yes" & length(config$cox_regression$multivariate_variables > 0)) {
  variables <- paste(config$cox_regression$multivariate_variables, collapse="+")
  multi_cox_formula <- as.formula(paste("Surv(time=patient_annotation[[config$patient_annotation_info$overall_survival_time_column_name]], 
                                 event=patient_annotation[[config$patient_annotation_info$overall_survival_event_column_name]])~", variables))
  multi_cox_fit <- coxph(multi_cox_formula, data=patient_annotation)
  print(multi_cox_fit)
  ggforest(multi_cox_fit)
} else if (config$cox_regression$do_overall_sex == "no") {
  print("You do not want a multivariate Cox regression analysis, so skipping this")
}


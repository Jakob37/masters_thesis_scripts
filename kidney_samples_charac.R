## FURTHER CHARACTERIZE KIDNEY SAMPLES -----------

library(tidyverse)
library(survival)
library(survminer)

kidney_class <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/data/kidney_clams_taxonomy.csv")
kidney_class <- kidney_class %>% mutate(cancer.type = substr(sample_id, 1, 4))
kidney_class$cancer.type <- factor(kidney_class$cancer.type)
kidney_class$clams.class <- factor(kidney_class$clams.class)
kidney_class$taxonomy_published <- factor(kidney_class$taxonomy_published)
kidney_taxonomy <- read_csv("/media/deboraholi/Data/LUND/9 THESIS/data/kidney_taxonomy.csv")
kidney_class <- left_join(kidney_class, kidney_taxonomy, by=c("taxonomy_published"="predom_enriched"))

kidney_class$taxonomy_published <- factor(kidney_class$taxonomy_published)
kidney_class$predom <- factor(kidney_class$predom)

rm(kidney_taxonomy)
# kidney_class_summary <- kidney_class %>% count(cancer.type, clams.class, predom, taxonomy_published)
# rm(kidney_class_summary)

addmargins(table(kidney_class$clams.class, kidney_class$predom, kidney_class$cancer.type, useNA = "ifany"))

# KICH
#         CC  Ch mixed   P <NA> Sum
# NonTRU   2  42     3   0    0  47
# TRU      0  14     1   0    2  17
# Sum      2  56     4   0    2  64

# KIRC
#         CC  Ch mixed   P <NA> Sum
# NonTRU 352  10     3   6    6 377
# TRU    114   5     7   7    2 135
# Sum    466  15    10  13    8 512

# KIRP
#         CC  Ch mixed   P <NA> Sum
# NonTRU   4   0     3 161   10 178
# TRU      1   2     4  95    3 105
# Sum      5   2     7 256   13 283

addmargins(table(kidney_class$clams.class, kidney_class$taxonomy_published, kidney_class$cancer.type, useNA = "ifany"))

# KICH
#        CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      0      0      2   42     3      0      0     0        0    0  47
# TRU         0      0      0   14     1      0      0     0        0    2  17
# Sum         0      0      2   56     4      0      0     0        0    2  64
# 
# KIRC
#        CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU     78    155    119   10     3      0      0     4        2    6 377
# TRU        19     91      4    5     7      2      1     4        0    2 135
# Sum        97    246    123   15    10      2      1     8        2    8 512
# 
# KIRP
#        CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      1      1      2    0     3     83     29    27       22   10 178
# TRU         0      1      0    2     4     43     39    13        0    3 105
# Sum         1      2      2    2     7    126     68    40       22   13 283


# REDIVIDE THEM ACCORDING TO DAVID'S CLASSIFICATION (FROM CHEN ET AL 2016) AND DO OS AND COX  -----------
addmargins(table(kidney_class$clams.class, kidney_class$predom, useNA = "ifany"))
#         CC  Ch mixed   P <NA> Sum
# NonTRU 358  52     9 167   16 602
# TRU    115  21    12 102    7 257
# Sum    473  73    21 269   23 859

addmargins(table(kidney_class$clams.class, kidney_class$taxonomy_published, kidney_class$predom, useNA = "ifany"))
# CC
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU     79    156    123    0     0      0      0     0        0    0 358
# TRU        19     92      4    0     0      0      0     0        0    0 115
# Sum        98    248    127    0     0      0      0     0        0    0 473
# 
# Ch
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      0      0      0   52     0      0      0     0        0    0  52
# TRU         0      0      0   21     0      0      0     0        0    0  21
# Sum         0      0      0   73     0      0      0     0        0    0  73
# 
# mixed
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      0      0      0    0     9      0      0     0        0    0   9
# TRU         0      0      0    0    12      0      0     0        0    0  12
# Sum         0      0      0    0    21      0      0     0        0    0  21
# 
# P
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      0      0      0    0     0     83     29    31       24    0 167
# TRU         0      0      0    0     0     45     40    17        0    0 102
# Sum         0      0      0    0     0    128     69    48       24    0 269
# 
# NA
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU      0      0      0    0     0      0      0     0        0   16  16
# TRU         0      0      0    0     0      0      0     0        0    7   7
# Sum         0      0      0    0     0      0      0     0        0   23  23
# 
# Sum
#         CC-e.1 CC-e.2 CC-e.3 Ch-e mixed P-e.1a P-e.1b P-e.2 P.CIMP-e <NA> Sum
# NonTRU     79    156    123   52     9     83     29    31       24   16 602
# TRU        19     92      4   21    12     45     40    17        0    7 257
# Sum        98    248    127   73    21    128     69    48       24   23 859



# OS and Cox
# --> get patient_annotation, censor, etc from OS.R first
patients_os_clams_kidney <- inner_join(patient_os_all_datasets, kidney_class,
                               by = c("sample.id" = "sample_id"))
patients_os_clams_kidney <- left_join(patients_os_clams_kidney, patient_annotation_prolif)

time_column <- "OS.time.years.5y"
event_column <- "OS.event.5y"
time_type <- "years"

groups.to.analyze.k <- levels(patients_os_clams_kidney$predom)

# CLAMS
# sink("kidney_david_os.txt") # to save output to a file
for (type_to_analyze in groups.to.analyze.k) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams_kidney, predom == type_to_analyze)
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
# sink()

# sink("kidney_david_cox.txt") # to save output to a file
for (type_to_analyze in groups.to.analyze.k) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams_kidney, predom == type_to_analyze)
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

# CLAMS by enrichment 
for (type_to_analyze in groups.to.analyze.k) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams_kidney, predom == type_to_analyze)
  type_subset$clams.class <- factor(type_subset$clams.class, levels=c("TRU", "NonTRU"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  # separating estimates by TRU or NonTRU
  fit_clams <- survfit(overall_surv_object~clams.class+taxonomy_published, data=type_subset)
  # testing by clams class
  print(survdiff(Surv(type_subset[[time_column]], 
                      type_subset[[event_column]])~type_subset$clams.class+type_subset$taxonomy_published,))
  
  # plot
  # simple for comparing info and making sure names are correct
  print(ggsurvplot(fit_clams, data=type_subset, title="Overall Survival", risk.table=TRUE))
  # good looking one
  # current_plot <- ggsurvplot(fit_clams, data=type_subset,
  #                            palette=c("darkorange1", "deepskyblue3"),
  #                            title=paste0("Overall Survival (", type_to_analyze, ")"), 
  #                            xlab=paste0("Time (", time_type, ")"),
  #                            censor.shape=124, censor.size=3,
  #                            pval=TRUE, pval.coord=c(0,0.1),
  #                            surv.median.line="hv",
  #                            risk.table=TRUE,
  #                            risk.table.fontsize = 4,
  #                            tables.theme = theme_survminer(font.main = 14),
  #                            legend="none", legend.title="CLAMS",
  #                            legend.labs=c("TRU", "NonTRU"))
  # print(current_plot)
  # current_filename <- paste0("OS 5y ", type_to_analyze, ".png")
  # ggsave(file=current_filename, print(current_plot), width=6.68, height=6.1, dpi=300)
}


# 25:75 PROLIFERATION
# redo Low/High in new cancer.type
patients_os_clams_kidney$prolif.kidney <- NA
for (cancer.group in groups.to.analyze.k) {
  print(cancer.group)
  # get only the data for that dataset + cancer type
  current.data <- subset(patients_os_clams_kidney, predom == cancer.group)
  # if karl.value < 25th percentile, call it "Low", else "High"
  current.data$prolif.group.karl.kidney <- ifelse(current.data$karl.value < quantile(current.data$karl.value, .25), "Low", "High")
  lows <- current.data[ which(current.data$prolif.group.karl.kidney == "Low"), ] %>% pull(sample.id)
  highs <- current.data[ which(current.data$prolif.group.karl.kidney == "High"), ] %>% pull(sample.id)
  patients_os_clams_kidney$prolif.kidney <- ifelse(patients_os_clams_kidney$sample.id %in% lows, "Low",
                                                        ifelse(patients_os_clams_kidney$sample.id %in% highs, "High", 
                                                               patients_os_clams_kidney$prolif.kidney))
}

# cox
for (type_to_analyze in groups.to.analyze.k) {
  print(type_to_analyze)
  type_subset <- subset(patients_os_clams_kidney, predom == type_to_analyze)
  type_subset$prolif.kidney <- factor(type_subset$prolif.kidney, levels=c("Low", "High"))
  # make the object
  overall_surv_object <- Surv(time=type_subset[[time_column]], 
                              event=type_subset[[event_column]])
  print(summary(coxph(formula = overall_surv_object~prolif.kidney, data = type_subset)))
  forest <- ggforest(coxph(formula = overall_surv_object~prolif.kidney, data = type_subset), 
                     main = paste("Hazard ratio", type_to_analyze))
  print(forest)
  current_filename <- paste0("Prolif Cox ", type_to_analyze, ".png")
  ggsave(file=current_filename, print(forest), width=6, height=3, dpi=300)
}

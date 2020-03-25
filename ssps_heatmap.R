#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO CREATE MATRIX FOR GENE RULES IN SSPS AND HEATMAP
#
#


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(pheatmap)
library(viridisLite)


# INPUT DATA FROM data_input.R  ----------------------------------

# (3) GEX
# filtered gex_matrix (only genes in ssp) from specific analyses gex_pam50.RData
# patient_annotation_ - classification results

# Functions  -------

create_empty_matrix <- function(original_gex_matrix, list_of_list_rules) {
  # create empty matrix for heatmap 
  # with rules as columns (extracted from aims.gs) 
  # and samples as rows (extracted from the original gex matrix)
  how_many_rules <- length(list_of_list_rules) * length(list_of_list_rules[[1]])
  rules_matrix <- matrix(nrow=length(colnames(original_gex_matrix)),
                         ncol=how_many_rules)
  row.names(rules_matrix) <- colnames(original_gex_matrix)
  rules_names <- c()
  for (rule in list_of_list_rules) {
    rules_names <- c(rules_names, rule)
  }
  colnames(rules_matrix) <- rules_names
  return(rules_matrix)
}

populate_matrix <- function(empty_matrix, original_gex_matrix, list_of_list_rules) {
  # now populate with increasing values for each list
  # if gene not in original gex matrix, put -1
  populated_matrix <- empty_matrix
  for (sample in row.names(empty_matrix)) {
    for (i in 1:length(list_of_list_rules)) {
      rules_list <- list_of_list_rules[[i]]
      for (rule in rules_list) {
        genes <- strsplit(rule, "<")
        genes <- unlist(genes)
        if (genes[1] %in% row.names(original_gex_matrix) & genes[2] %in% row.names(original_gex_matrix)) {
          ifelse(original_gex_matrix[genes[1],sample] < original_gex_matrix[genes[2],sample],
                 populated_matrix[sample,rule] <- i,
                 populated_matrix[sample,rule] <- 0)
        } else {
          populated_matrix[sample,rule] <- -1
        }
      }
    }
  }
  return(populated_matrix)
}



# CLAMS  ---------------------------------------------

library(CLAMS)

clams_rules_tru <- CLAMSmodel$selected.pairs.list$bronchioid
clams_rules_nontru <- CLAMSmodel$selected.pairs.list$nonbronchioid

clams_rules <- list(clams_rules_tru=clams_rules_tru, clams_rules_nontru=clams_rules_nontru)

annotation_colors = list(
  clams.class = c(TRU = "#FF7F00", NonTRU = "#009ACD")
)


# GOBO

gobo_empty_matrix <- create_empty_matrix(gex_matrix_gobo_clams, clams_rules)
gobo_matrix_to_plot <- populate_matrix(gobo_empty_matrix, gex_matrix_gobo_clams, clams_rules)

row_annotation_gobo <- patient_annotation_clams %>% filter(dataset == "GOBO")
row_annotation_gobo <- row_annotation_gobo[c("sample.id", "clams.class")]
row_annotation_gobo$clams.class <- factor(row_annotation_gobo$clams.class, levels=c("TRU", "NonTRU"))
row_annotation_gobo <- row_annotation_gobo %>% arrange(clams.class)
row_annotation_gobo <- column_to_rownames(row_annotation_gobo, var = "sample.id")

clams_gobo_matrix <- clams_gobo_matrix[rownames(row_annotation_gobo), ]

gobo_heatmap <- pheatmap(mat = clams_gobo_matrix, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_gobo,
                         #annotation_col = col_annotation,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         color = c("white", "darkorange1", "deepskyblue3"),
                         main = "GOBO, CLAMS"
)


# SCAN-B

scanb_empty_matrix <- create_empty_matrix(gex_matrix_scanb_clams, clams_rules)
scanb_matrix_to_plot <- populate_matrix(scanb_empty_matrix, gex_matrix_scanb_clams, clams_rules)

# to be done if needed

# ALL TCGA

tcga_empty_matrix <- create_empty_matrix(gex_matrix_tcga_clams, clams_rules)
tcga_matrix_to_plot <- populate_matrix(tcga_empty_matrix, gex_matrix_tcga_clams, clams_rules)

# to be done if needed

# ALL SAMPLES

all_matrix_to_plot <- rbind(gobo_matrix_to_plot, scanb_matrix_to_plot)
all_matrix_to_plot <- rbind(all_matrix_to_plot, tcga_matrix_to_plot)

# create class annotation
row_annotation_tcga <- patient_annotation_clams[c("sample.id", "clams.class")]
row_annotation_tcga$clams.class <- factor(row_annotation_tcga$clams.class, c("TRU", "NonTRU"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(clams.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

all_matrix_to_plot <- all_matrix_to_plot[rownames(row_annotation_tcga), ]

all_tcga_heatmap <- pheatmap(mat = all_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation_tcga,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              annotation_names_col = FALSE,
                              color = c("white", "darkorange1", "deepskyblue3"),
                              main = "All datasets, CLAMS"
)


# TCGA LUAD

luad_tcga_subset_samples <- subset(patient_annotation_clams, dataset == 'TCGA' & cancer.type == "LUAD") %>% pull(sample.id)
luad_tcga_gex_matrix <- gex_matrix_tcga_clams[,luad_tcga_subset_samples]

luad_tcga_empty_matrix <- create_empty_matrix(luad_tcga_gex_matrix, clams_rules)
luad_tcga_matrix_to_plot <- populate_matrix(luad_tcga_empty_matrix, luad_tcga_gex_matrix, clams_rules)

# create class annotation
row_annotation_tcga <- patient_annotation_clams %>% filter(dataset == "TCGA" & cancer.type == "LUAD")
row_annotation_tcga <- row_annotation_tcga[c("sample.id", "clams.class")]
row_annotation_tcga$clams.class <- factor(row_annotation_tcga$clams.class, c("TRU", "NonTRU"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(clams.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

luad_tcga_matrix_to_plot <- luad_tcga_matrix_to_plot[rownames(row_annotation_tcga), ]

luad_tcga_heatmap <- pheatmap(mat = luad_tcga_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation_tcga,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              annotation_names_col = FALSE,
                              color = c("white", "darkorange1", "deepskyblue3"),
                              main = "TCGA LUAD, CLAMS"
)

# only TRU rules
luad_tcga_empty_matrix <- create_empty_matrix(luad_tcga_gex_matrix, clams_rules_tru)
luad_tcga_matrix_to_plot <- populate_matrix(luad_tcga_empty_matrix, luad_tcga_gex_matrix, clams_rules_tru)
luad_tcga_matrix_to_plot[luad_tcga_matrix_to_plot > 1] <- 1

# create class annotation
row_annotation_tcga <- patient_annotation_clams %>% filter(dataset == "TCGA" & cancer.type == "LUAD")
row_annotation_tcga <- row_annotation_tcga[c("sample.id", "clams.class")]
row_annotation_tcga$clams.class <- factor(row_annotation_tcga$clams.class, c("TRU", "NonTRU"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(clams.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

luad_tcga_matrix_to_plot <- luad_tcga_matrix_to_plot[rownames(row_annotation_tcga), ]

luad_tcga_heatmap <- pheatmap(mat = luad_tcga_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                              show_colnames = TRUE,
                              gaps_row = c(152),
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation_tcga,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              annotation_names_col = FALSE,
                              color = c("white", "grey70"),
                              main = "TCGA LUAD, CLAMS"
)

# ALL BRCA together

brca_tcga_subset_samples <- subset(patient_annotation_clams, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
brca_tcga_gex_matrix <- gex_matrix_tcga_clams[,brca_tcga_subset_samples]

brca_tcga_empty_matrix <- create_empty_matrix(brca_tcga_gex_matrix, clams_rules_tru)
brca_tcga_matrix_to_plot <- populate_matrix(brca_tcga_empty_matrix, brca_tcga_gex_matrix, clams_rules_tru)
brca_tcga_matrix_to_plot[brca_tcga_matrix_to_plot > 1] <- 1

gobo_empty_matrix <- create_empty_matrix(gex_matrix_gobo_clams, clams_rules_tru)
gobo_matrix_to_plot <- populate_matrix(gobo_empty_matrix, gex_matrix_gobo_clams, clams_rules_tru)
gobo_matrix_to_plot[gobo_matrix_to_plot > 1] <- 1

scanb_empty_matrix <- create_empty_matrix(gex_matrix_scanb_clams, clams_rules_tru)
scanb_matrix_to_plot <- populate_matrix(scanb_empty_matrix, gex_matrix_scanb_clams, clams_rules_tru)
scanb_matrix_to_plot[scanb_matrix_to_plot > 1] <- 1

all_brca_matrix_to_plot <- rbind(brca_tcga_matrix_to_plot, gobo_matrix_to_plot, scanb_matrix_to_plot)

# create class annotation
row_annotation <- patient_annotation_clams %>% filter(cancer.type == "BRCA")
row_annotation <- row_annotation[c("sample.id", "clams.class", "dataset")]
row_annotation$clams.class <- factor(row_annotation$clams.class, c("TRU", "NonTRU"))
row_annotation <- row_annotation %>% arrange(clams.class)
row_annotation$dataset <- factor(row_annotation$dataset, c("TCGA", "GOBO", "SCAN-B"))
row_annotation <- row_annotation %>% arrange(dataset)
row_annotation <- row_annotation[-3]
row_annotation <- column_to_rownames(row_annotation, var = "sample.id")

all_brca_matrix_to_plot <- all_brca_matrix_to_plot[rownames(row_annotation), ]

all_brca_heatmap <- pheatmap(mat = all_brca_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                             show_colnames = TRUE,
                              gaps_row = c(1072,2953),
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              color = c("white", "grey70"),
                              main = "ALL BRCA, CLAMS"
)


# PAM50  --------------------------------------------------

# Reduced

pam50_rules_luma <- pam50.red.aims.gs$selected.pairs.list$LumA
pam50_rules_lumb <- pam50.red.aims.gs$selected.pairs.list$LumB
pam50_rules_her2 <- pam50.red.aims.gs$selected.pairs.list$Her2
pam50_rules_basal <- pam50.red.aims.gs$selected.pairs.list$Basal

pam50_red_rules <- list(pam50_rules_luma=pam50_rules_luma, pam50_rules_lumb=pam50_rules_lumb, 
                    pam50_rules_her2=pam50_rules_her2, pam50_rules_basal=pam50_rules_basal)

annotation_colors = list(
  pam50.red.class = c(LumA = "#F0F921FF", LumB = "#F99B3EFF", Her2 = "#D35171FF", Basal = "#900DA4FF")
)


# GOBO

gobo_empty_matrix <- create_empty_matrix(gex_matrix_gobo_red, pam50_rules)
gobo_matrix_to_plot <- populate_matrix(gobo_empty_matrix, gex_matrix_gobo_red, pam50_rules)

# create class annotation
row_annotation_gobo <- patient_annotation_pam50 %>% filter(dataset == "GOBO")
row_annotation_gobo <- row_annotation_gobo[c("sample.id", "pam50.red.class")]
row_annotation_gobo$pam50.red.class <- factor(row_annotation_gobo$pam50.red.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_gobo <- row_annotation_gobo %>% arrange(pam50.red.class)
row_annotation_gobo <- column_to_rownames(row_annotation_gobo, var = "sample.id")

gobo_matrix_to_plot <- gobo_matrix_to_plot[rownames(row_annotation_gobo), ]

gobo_heatmap <- pheatmap(mat = gobo_matrix_to_plot, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         gaps_row = c(615,867,1656),
                         gaps_col = c(22,44,66),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_gobo,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                         main = "GOBO, PAM50, reduced"
)


# SCAN-B

scanb_empty_matrix <- create_empty_matrix(gex_matrix_scanb_red, pam50_rules)
scanb_matrix_to_plot <- populate_matrix(scanb_empty_matrix, gex_matrix_scanb_red, pam50_rules)

# create class annotation
row_annotation_scanb <- patient_annotation_pam50 %>% filter(dataset == "SCAN-B")
row_annotation_scanb <- row_annotation_scanb[c("sample.id", "pam50.red.class")]
row_annotation_scanb$pam50.red.class <- factor(row_annotation_scanb$pam50.red.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_scanb <- row_annotation_scanb %>% arrange(pam50.red.class)
row_annotation_scanb <- column_to_rownames(row_annotation_scanb, var = "sample.id")

scanb_matrix_to_plot <- scanb_matrix_to_plot[rownames(row_annotation_scanb), ]

the_heatmap <- pheatmap(mat = scanb_matrix_to_plot, 
                        border_color = NA, 
                        show_rownames = FALSE,
                        show_colnames = FALSE,
                        gaps_row = c(1842,2693,3165),
                        gaps_col = c(22,44,66),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        annotation_row = row_annotation_scanb,
                        annotation_colors = annotation_colors,
                        annotation_names_row = FALSE,
                        annotation_names_col = FALSE,
                        color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                        main = "SCAN-B, PAM50, reduced"
)


# BRCA TCGA

brca_tcga_subset_samples <- subset(patient_annotation_pam50, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
brca_tcga_gex_matrix <- gex_matrix_tcga_red[,brca_tcga_subset_samples]

brca_tcga_empty_matrix <- create_empty_matrix(brca_tcga_gex_matrix, pam50_rules)
brca_tcga_matrix_to_plot <- populate_matrix(brca_tcga_empty_matrix, brca_tcga_gex_matrix, pam50_rules)

# create class annotation
row_annotation_tcga <- patient_annotation_pam50 %>% filter(dataset == "TCGA" & cancer.type == "BRCA")
row_annotation_tcga <- row_annotation_tcga[c("sample.id", "pam50.red.class")]
row_annotation_tcga$pam50.red.class <- factor(row_annotation_tcga$pam50.red.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(pam50.red.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

brca_tcga_matrix_to_plot <- brca_tcga_matrix_to_plot[rownames(row_annotation_tcga), ]

brca_tcga_heatmap <- pheatmap(mat = brca_tcga_matrix_to_plot, 
                               border_color = NA, 
                               show_rownames = FALSE,
                               show_colnames = FALSE,
                               gaps_row = c(255,661,885),
                               gaps_col = c(22,44,66),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               annotation_row = row_annotation_tcga,
                               annotation_colors = annotation_colors,
                               annotation_names_row = FALSE,
                               annotation_names_col = FALSE,
                               color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                               main = "TCGA BRCA, PAM50, reduced"
)


# Full set

pam50_rules_luma <- pam50.all.aims.gs$selected.pairs.list$LumA
pam50_rules_lumb <- pam50.all.aims.gs$selected.pairs.list$LumB
pam50_rules_her2 <- pam50.all.aims.gs$selected.pairs.list$Her2
pam50_rules_basal <- pam50.all.aims.gs$selected.pairs.list$Basal

pam50_all_rules <- list(pam50_rules_luma=pam50_rules_luma, pam50_rules_lumb=pam50_rules_lumb, 
                        pam50_rules_her2=pam50_rules_her2, pam50_rules_basal=pam50_rules_basal)

annotation_colors = list(
  pam50.all.class = c(LumA = "#F0F921FF", LumB = "#F99B3EFF", Her2 = "#D35171FF", Basal = "#900DA4FF")
)


# GOBO

gobo_empty_matrix <- create_empty_matrix(gex_matrix_gobo_pam50, pam50_all_rules)
gobo_matrix_to_plot <- populate_matrix(gobo_empty_matrix, gex_matrix_gobo_pam50, pam50_all_rules)

# create class annotation
row_annotation_gobo <- patient_annotation_pam50 %>% filter(dataset == "GOBO")
row_annotation_gobo <- row_annotation_gobo[c("sample.id", "pam50.all.class")]
row_annotation_gobo$pam50.all.class <- factor(row_annotation_gobo$pam50.all.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_gobo <- row_annotation_gobo %>% arrange(pam50.all.class)
row_annotation_gobo <- column_to_rownames(row_annotation_gobo, var = "sample.id")

gobo_matrix_to_plot <- gobo_matrix_to_plot[rownames(row_annotation_gobo), ]

gobo_heatmap <- pheatmap(mat = gobo_matrix_to_plot, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         gaps_row = c(505,1165,1610),
                         gaps_col = c(21,42,63),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_gobo,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         color = c("grey80", "white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                         main = "GOBO, PAM50, full set"
)

# removing all the missing rules
gobo_matrix_to_plot_clean <- as_tibble(gobo_matrix_to_plot, rownames=NA) %>% select_if(~sum(. > 0) > 0)

gobo_heatmap_clean <- pheatmap(mat = gobo_matrix_to_plot_clean, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         gaps_row = c(505,1165,1610),
                         gaps_col = c(9,21,35),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_gobo,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                         main = "GOBO, PAM50, full set, without missing rules"
)

# SCAN-B

scanb_empty_matrix <- create_empty_matrix(gex_matrix_scanb_pam50, pam50_all_rules)
scanb_matrix_to_plot <- populate_matrix(scanb_empty_matrix, gex_matrix_scanb_pam50, pam50_all_rules)

# create class annotation
row_annotation_scanb <- patient_annotation_pam50 %>% filter(dataset == "SCAN-B")
row_annotation_scanb <- row_annotation_scanb[c("sample.id", "pam50.all.class")]
row_annotation_scanb$pam50.all.class <- factor(row_annotation_scanb$pam50.all.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_scanb <- row_annotation_scanb %>% arrange(pam50.all.class)
row_annotation_scanb <- column_to_rownames(row_annotation_scanb, var = "sample.id")

scanb_matrix_to_plot <- scanb_matrix_to_plot[rownames(row_annotation_scanb), ]

scanb_heatmap <- pheatmap(mat = scanb_matrix_to_plot, 
                        border_color = NA, 
                        show_rownames = FALSE,
                        show_colnames = FALSE,
                        gaps_row = c(1822,2636,3150),
                        gaps_col = c(21,42,63),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        annotation_row = row_annotation_scanb,
                        annotation_colors = annotation_colors,
                        annotation_names_row = FALSE,
                        annotation_names_col = FALSE,
                        color = c("grey80", "white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                        main = "SCAN-B, PAM50, full set"
)

# removing all the missing rules
scanb_matrix_to_plot_clean <- as_tibble(scanb_matrix_to_plot, rownames=NA) %>% select_if(~sum(. > 0) > 0)

scanb_heatmap_clean <- pheatmap(mat = scanb_matrix_to_plot_clean, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         gaps_row = c(1822,2636,3150),
                         gaps_col = c(21,41,61),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_scanb,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                         main = "SCAN-B, PAM50, full set, without missing rules"
)



# BRCA TCGA

brca_tcga_subset_samples <- subset(patient_annotation_pam50, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
brca_tcga_gex_matrix <- gex_matrix_tcga_pam50[,brca_tcga_subset_samples]

brca_tcga_empty_matrix <- create_empty_matrix(brca_tcga_gex_matrix, pam50_all_rules)
brca_tcga_matrix_to_plot <- populate_matrix(brca_tcga_empty_matrix, brca_tcga_gex_matrix, pam50_all_rules)

# create class annotation
row_annotation_tcga <- patient_annotation_pam50 %>% filter(dataset == "TCGA" & cancer.type == "BRCA")
row_annotation_tcga <- row_annotation_tcga[c("sample.id", "pam50.all.class")]
row_annotation_tcga$pam50.all.class <- factor(row_annotation_tcga$pam50.all.class, c("LumA", "LumB", "Her2", "Basal"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(pam50.all.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

brca_tcga_matrix_to_plot <- brca_tcga_matrix_to_plot[rownames(row_annotation_tcga), ]

brca_tcga_heatmap <- pheatmap(mat = brca_tcga_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              gaps_row = c(251,643,880),
                              gaps_col = c(21,42,63),
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation_tcga,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              annotation_names_col = FALSE,
                              color = c("white", "#F0F921FF", "#F99B3EFF", "#D35171FF", "#900DA4FF"),
                              main = "TCGA BRCA, PAM50, full set"
)


# ROR  --------------------------------------------------


annotation_colors = list(
  ror.red.hl.class = c(Low = "#FDE725FF", Medium = "#23898EFF", High = "#440154FF"),
  ror.red.class = c(c005="#FDE725FF", c010="#D8E219FF", c015="#B0DD2FFF", c020="#89D548FF",
                c025="#65CB5EFF", c030="#46C06FFF", c035="#2EB37CFF", c040="#21A585FF",
                c045="#1F978BFF", c050="#23898EFF", c055="#297B8EFF", c060="#2E6D8EFF",
                c065="#355E8DFF", c070="#3D4E8AFF", c075="#433D84FF", c080="#472A7AFF",
                c085="#481769FF", c090="#440154FF") # from viridis(18), descending
)



# Reduced

ror_rules_c005 <- ror.red.aims.gs$selected.pairs.list$c005
ror_rules_c010 <- ror.red.aims.gs$selected.pairs.list$c010
ror_rules_c015 <- ror.red.aims.gs$selected.pairs.list$c015
ror_rules_c020 <- ror.red.aims.gs$selected.pairs.list$c020
ror_rules_c025 <- ror.red.aims.gs$selected.pairs.list$c025
ror_rules_c030 <- ror.red.aims.gs$selected.pairs.list$c030
ror_rules_c035 <- ror.red.aims.gs$selected.pairs.list$c035
ror_rules_c040 <- ror.red.aims.gs$selected.pairs.list$c040
ror_rules_c045 <- ror.red.aims.gs$selected.pairs.list$c045
ror_rules_c050 <- ror.red.aims.gs$selected.pairs.list$c050
ror_rules_c055 <- ror.red.aims.gs$selected.pairs.list$c055
ror_rules_c060 <- ror.red.aims.gs$selected.pairs.list$c060
ror_rules_c065 <- ror.red.aims.gs$selected.pairs.list$c065
ror_rules_c070 <- ror.red.aims.gs$selected.pairs.list$c070
ror_rules_c075 <- ror.red.aims.gs$selected.pairs.list$c075
ror_rules_c080 <- ror.red.aims.gs$selected.pairs.list$c080
ror_rules_c085 <- ror.red.aims.gs$selected.pairs.list$c085
ror_rules_c090 <- ror.red.aims.gs$selected.pairs.list$c090

ror_red_rules <- list(
  ror_rules_c005 = ror_rules_c005,
  ror_rules_c010 = ror_rules_c010,
  ror_rules_c015 = ror_rules_c015,
  ror_rules_c020 = ror_rules_c020,
  ror_rules_c025 = ror_rules_c025,
  ror_rules_c030 = ror_rules_c030,
  ror_rules_c035 = ror_rules_c035,
  ror_rules_c040 = ror_rules_c040,
  ror_rules_c045 = ror_rules_c045,
  ror_rules_c050 = ror_rules_c050,
  ror_rules_c055 = ror_rules_c055,
  ror_rules_c060 = ror_rules_c060,
  ror_rules_c065 = ror_rules_c065,
  ror_rules_c070 = ror_rules_c070,
  ror_rules_c075 = ror_rules_c075,
  ror_rules_c080 = ror_rules_c080,
  ror_rules_c085 = ror_rules_c085,
  ror_rules_c090 = ror_rules_c090
)



# GOBO

gobo_empty_matrix <- create_empty_matrix(gex_matrix_gobo_ror, ror_red_rules)
gobo_matrix_to_plot <- populate_matrix(gobo_empty_matrix, gex_matrix_gobo_ror, ror_red_rules)

# create class annotation
row_annotation_gobo <- patient_annotation_ror %>% filter(dataset == "GOBO")
row_annotation_gobo <- row_annotation_gobo[c("sample.id", "ror.red.class", "ror.red.hl.class")]
row_annotation_gobo$ror.red.class <- factor(row_annotation_gobo$ror.red.class, 
                                            levels = c("c005", "c010", "c015", "c020", "c025", "c030",
                                                      "c035", "c040", "c045", "c050", "c055",
                                                      "c060", "c065", "c070", "c075", "c080",
                                                      "c085", "c090"))
row_annotation_gobo <- row_annotation_gobo %>% arrange(ror.red.class)
row_annotation_gobo <- column_to_rownames(row_annotation_gobo, var = "sample.id")

gobo_matrix_to_plot <- gobo_matrix_to_plot[rownames(row_annotation_gobo), ]

# get colors of only classes that are there
row_annotation_gobo %>% group_by(ror.red.class) %>% tally()

gobo_heatmap <- pheatmap(mat = gobo_matrix_to_plot, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         #gaps_row = c(437,524,884,1091,1502,1617,1655,1665,1666),
                         #gaps_col = c(),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_gobo,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         color = c("white", "#FDE725FF", "#21A585FF", "#1F978BFF", "#23898EFF",
                                   "#297B8EFF", "#2E6D8EFF", "#355E8DFF", "#3D4E8AFF",
                                   "#481769FF", "#440154FF"),
                         main = "GOBO, ROR, reduced"
)


# SCAN-B

scanb_empty_matrix <- create_empty_matrix(gex_matrix_scanb_ror, ror_red_rules)
scanb_matrix_to_plot <- populate_matrix(scanb_empty_matrix, gex_matrix_scanb_ror, ror_red_rules)

# create class annotation
row_annotation_scanb <- patient_annotation_ror %>% filter(dataset == "SCAN-B")
row_annotation_scanb <- row_annotation_scanb[c("sample.id", "ror.red.class", "ror.red.hl.class")]
row_annotation_scanb$ror.red.class <- factor(row_annotation_scanb$ror.red.class, 
                                            levels = c("c005", "c010", "c015", "c020", "c025", "c030",
                                                       "c035", "c040", "c045", "c050", "c055",
                                                       "c060", "c065", "c070", "c075", "c080",
                                                       "c085", "c090"))
row_annotation_scanb <- row_annotation_scanb %>% arrange(ror.red.class)
row_annotation_scanb <- column_to_rownames(row_annotation_scanb, var = "sample.id")

scanb_matrix_to_plot <- scanb_matrix_to_plot[rownames(row_annotation_scanb), ]

# get colors of only classes that are there
row_annotation_scanb %>% group_by(ror.red.class) %>% tally()

scanb_heatmap <- pheatmap(mat = scanb_matrix_to_plot, 
                         border_color = NA, 
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         #gaps_row = c(1213,1215,1230,1237,1421,1567,
                          #            1669,1873,2042,2298,2516,2662,
                          #            2849,2926,2973,2978),
                         #gaps_col = c(),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = row_annotation_scanb,
                         annotation_colors = annotation_colors,
                         annotation_names_row = FALSE,
                         color = c("grey80", "white", "white", "#FDE725FF", "#B0DD2FFF", "#89D548FF",
                                   "#65CB5EFF", "#46C06FFF", "#2EB37CFF", "#21A585FF",
                                   "#1F978BFF", "#23898EFF", "#297B8EFF", "#2E6D8EFF",
                                   "#355E8DFF", "#3D4E8AFF", "#433D84FF", "#472A7AFF",
                                   "#481769FF", "#440154FF"),
                         main = "SCAN-B, ROR, reduced"
)

# removing all the missing rules
scanb_matrix_to_plot_clean <- as_tibble(scanb_matrix_to_plot, rownames=NA) %>% select_if(~sum(. > 0) > 0)

scanb_heatmap_clean <- pheatmap(mat = scanb_matrix_to_plot_clean, 
                          border_color = NA, 
                          show_rownames = FALSE,
                          show_colnames = FALSE,
                          #gaps_row = c(1213,1215,1230,1237,1421,1567,
                          #            1669,1873,2042,2298,2516,2662,
                          #            2849,2926,2973,2978),
                          #gaps_col = c(),
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          annotation_row = row_annotation_scanb,
                          annotation_colors = annotation_colors,
                          annotation_names_row = FALSE,
                          color = c("grey80", "white", "#FDE725FF", "#89D548FF",
                                    "#46C06FFF", "#2EB37CFF", "#21A585FF",
                                    "#1F978BFF", "#23898EFF", "#297B8EFF", "#2E6D8EFF",
                                    "#355E8DFF", "#3D4E8AFF", "#433D84FF", "#472A7AFF",
                                    "#481769FF", "#440154FF"),
                          main = "SCAN-B, ROR, reduced, no missing rules"
)

# BRCA TCGA

brca_tcga_subset_samples <- subset(patient_annotation_ror, dataset == 'TCGA' & cancer.type == "BRCA") %>% pull(sample.id)
brca_tcga_gex_matrix <- gex_matrix_tcga_ror[,brca_tcga_subset_samples]

brca_tcga_empty_matrix <- create_empty_matrix(brca_tcga_gex_matrix, ror_red_rules)
brca_tcga_matrix_to_plot <- populate_matrix(brca_tcga_empty_matrix, brca_tcga_gex_matrix, ror_red_rules)

# create class annotation
row_annotation_tcga <- patient_annotation_ror %>% filter(dataset == "TCGA" & cancer.type == "BRCA")
row_annotation_tcga <- row_annotation_tcga[c("sample.id", "ror.red.class", "ror.red.hl.class")]
row_annotation_tcga$ror.red.class <- factor(row_annotation_tcga$ror.red.class, 
                                            levels = c("c005", "c010", "c015", "c020", "c025", "c030",
                                                       "c035", "c040", "c045", "c050", "c055",
                                                       "c060", "c065", "c070", "c075", "c080",
                                                       "c085", "c090"))
row_annotation_tcga <- row_annotation_tcga %>% arrange(ror.red.class)
row_annotation_tcga <- column_to_rownames(row_annotation_tcga, var = "sample.id")

brca_tcga_matrix_to_plot <- brca_tcga_matrix_to_plot[rownames(row_annotation_tcga), ]

# get colors of only classes that are there
row_annotation_tcga %>% group_by(ror.red.class) %>% tally()

brca_tcga_heatmap <- pheatmap(mat = brca_tcga_matrix_to_plot, 
                              border_color = NA, 
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              #gaps_row = c(113,115,121,125,161,200,260,
                              #             351,417,471,516,540,544,585),
                              #gaps_col = c(),
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_row = row_annotation_tcga,
                              annotation_colors = annotation_colors,
                              annotation_names_row = FALSE,
                              color = c("white", "#FDE725FF", "#89D548FF",
                                        "#46C06FFF", "#2EB37CFF", "#21A585FF",
                                        "#1F978BFF", "#23898EFF", "#297B8EFF", "#2E6D8EFF",
                                        "#355E8DFF", "#3D4E8AFF", "#433D84FF", "#472A7AFF",
                                        "#481769FF", "#440154FF"),
                              main = "TCGA BRCA, ROR, reduced"
)


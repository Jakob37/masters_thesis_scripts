#
# Create a report for the 3 datasets I'll be working with
#
#


# LOAD GENERAL PACKAGES ------------------------------------------

library(tidyverse)
library(viridisLite)
library(org.Hs.eg.db)


# INPUT DATA FROM data_input.R  ----------------------------------

# (2) GENE TABLE INFORMATION
# (3) GEX
# functions.R


## GOBO  ----------------------------------------------------

# Reporters
reporter_ids_gobo <- gex_matrix_gobo[,1:2] 
reporter_ids_gobo <- rownames_to_column(data.frame(reporter_ids_gobo), var = "reporterId")
reporter_ids_gobo <- reporter_ids_gobo$reporterId

sum(is.na(reporter_ids_gobo)) # 0 NAs, but test:
length(reporter_ids_gobo[reporter_ids_gobo == "NA"]) # 0 NAs
sum(!is.na(reporter_ids_gobo)) # 18223 not NAs
length(unique(reporter_ids_gobo)) # 18223 unique reporters present
sum(duplicated(reporter_ids_gobo)) # no duplicates

# EntrezID
entrez_ids_gobo <- gex_matrix_gobo[,1:2] 
entrez_ids_gobo <- rownames_to_column(data.frame(entrez_ids_gobo), var = "reporterId")
entrez_ids_gobo <- left_join(entrez_ids_gobo, gene_table_gobo[,c("reporterId","entrezId")])
entrez_ids_gobo <- entrez_ids_gobo %>% mutate(entrezId = paste0("e",entrezId))
entrez_ids_gobo <- entrez_ids_gobo$entrezId

sum(is.na(entrez_ids_gobo)) # 0 NAs, but test:
length(entrez_ids_gobo[entrez_ids_gobo == "NA"]) # 0 NAs
sum(!is.na(entrez_ids_gobo)) # 18223 not NAs
length(unique(entrez_ids_gobo)) # 12201 genes present
sum(duplicated(entrez_ids_gobo)) # 6022 duplicates
duplicate_gobo <- unique(entrez_ids_gobo[duplicated(entrez_ids_gobo)])
length(duplicate_gobo) # 4030 genes that are duplicated
gobo_dups_count <- get.dup.genes.count(entrez_ids_gobo)
sum(gobo_dups_count) # 18223, sanity check
as.data.frame(gobo_dups_count) %>% group_by(gobo_dups_count) %>% tally() %>% print(n=Inf)
# gobo_dups_count     n
# 1                 8171
# 2                 2639
# 3                 975
# 4                  306
# 5                   76
# 6                   18
# 7                   10
# 8                  1
# 9                  3
# 10                  1
# 22                  1

ggplot(mapping = aes(gobo_dups_count)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  theme_minimal() +
  scale_x_continuous(breaks = c(1:10,22)) +
  labs(x = "How many times the entrezId is in the table", y ="EntrezId count percentage") +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", 
            hjust=-0.2, vjust = 0.5, size=3, angle=90)
# 8171 ids (67%) appear just once in the table. 6 ids appear more than 8 times.



## SCAN-B  ----------------------------------------------------

# Gene symbols
gene_symbols_scanb <- gex_matrix_scanb[,1:2]
gene_symbols_scanb <- rownames_to_column(data.frame(gene_symbols_scanb), var = "geneSymbol")
gene_symbols_scanb <- gene_symbols_scanb$geneSymbol

sum(is.na(gene_symbols_scanb)) # 0 NAs, but test:
length(gene_symbols_scanb[gene_symbols_scanb == "NA"]) # 0 NAs
sum(!is.na(gene_symbols_scanb)) # 19102 not NAs
length(unique(gene_symbols_scanb)) # 19102 unique genes present
sum(duplicated(gene_symbols_scanb)) # no duplicates

# EntrezID
entrez_ids_scanb <- mapIds(org.Hs.eg.db, gene_symbols_scanb, 'ENTREZID', 'SYMBOL')
entrez_ids_scanb <- unname(entrez_ids_scanb)

sum(is.na(entrez_ids_scanb)) # 1030 NAs, but test:
length(entrez_ids_scanb[entrez_ids_scanb == "NA"]) # 0 NAs
sum(!is.na(entrez_ids_scanb)) # 18072 not NAs
length(unique(entrez_ids_scanb)) # 18072 genes present (disregard 1 NA from 18073)
sum(duplicated(entrez_ids_scanb)) # 1029 duplicates (NAs)
unique(entrez_ids_scanb[duplicated(entrez_ids_scanb)]) # just NA



## TCGA  ----------------------------------------------------

# ensemblID
ensembl_ids_tcga <- gex_matrix_tcga[,1:2]
ensembl_ids_tcga <- rownames_to_column(data.frame(ensembl_ids_tcga), var = "ensemblId")
ensembl_ids_tcga <- ensembl_ids_tcga$ensemblId

sum(is.na(ensembl_ids_tcga)) # 0 NAs, but test:
length(ensembl_ids_tcga[ensembl_ids_tcga == "NA"]) # 0 NAs
sum(!is.na(ensembl_ids_tcga)) # 19676 not NAs
length(unique(ensembl_ids_tcga)) # 19676 unique ensembl ids present
sum(duplicated(ensembl_ids_tcga)) # no duplicates

# EntrezID using gene annotation table that came with it in folder
ensembl_ids_tcga <- gex_matrix_tcga[,1:2]
ensembl_ids_tcga <- rownames_to_column(data.frame(ensembl_ids_tcga), var = "ensemblId")
  entrez_ids_table_tcga <- left_join(ensembl_ids_tcga, gene_table_tcga[,c("ENSG","ENTREZID")], by=c("ensemblId" = "ENSG"))
entrez_ids_table_tcga <- entrez_ids_table_tcga$ENTREZID

sum(is.na(entrez_ids_table_tcga)) # 0 NAs, but test:
length(entrez_ids_table_tcga[entrez_ids_table_tcga == "NA"]) # 579 NAs
sum(!is.na(entrez_ids_table_tcga)) # 19676 not NAs, but test:
length(entrez_ids_table_tcga[!entrez_ids_table_tcga == "NA"]) # 19097 not NAs
length(unique(entrez_ids_table_tcga[!entrez_ids_table_tcga == "NA"])) # 19081 unique genes present
sum(duplicated(entrez_ids_table_tcga)) # 594 duplicates, 579 NAs
duplicate_tcga <- unique(entrez_ids_table_tcga[duplicated(entrez_ids_table_tcga)])
length(duplicate_tcga) # 15 genes that are duplicated, excluding one NA
table_tcga_dups_count <- get.dup.genes.count(entrez_ids_table_tcga[!entrez_ids_table_tcga == "NA"])
sum(table_tcga_dups_count) # 20000, sanity check - but why?? should be 19097... because of the ones that have | in it... fix.
as.data.frame(table_tcga_dups_count) %>% group_by(table_tcga_dups_count) %>% tally() %>% print(n=Inf)

# EntrezID using org.Hs.eg.db
ensembl_ids_tcga <- gex_matrix_tcga[,1:2]
ensembl_ids_tcga <- rownames_to_column(data.frame(ensembl_ids_tcga), var = "ensemblId")
ensembl_ids_tcga <- ensembl_ids_tcga$ensemblId
entrez_ids_org_tcga <- mapIds(org.Hs.eg.db, ensembl_ids_tcga, 'ENTREZID', 'ENSEMBL')
entrez_ids_org_tcga <- unname(entrez_ids_org_tcga)

sum(is.na(entrez_ids_org_tcga)) # 453 NAs, but test:
length(entrez_ids_org_tcga[entrez_ids_org_tcga == "NA"]) # 453 NAs
sum(!is.na(entrez_ids_org_tcga)) # 19223 not NAs
length(unique(entrez_ids_org_tcga)) # 19177 unique genes present
sum(duplicated(entrez_ids_org_tcga)) # 499 duplicates
duplicate_tcga <- unique(entrez_ids_org_tcga[duplicated(entrez_ids_org_tcga)])
length(duplicate_tcga) # 44 genes that are duplicated, excluding one NA
org_tcga_dups_count <- get.dup.genes.count(entrez_ids_org_tcga[!is.na(entrez_ids_org_tcga)])
sum(org_tcga_dups_count) # 19223, sanity check
as.data.frame(org_tcga_dups_count) %>% group_by(org_tcga_dups_count) %>% tally() %>% print(n=Inf)
# 19132 ids appear just once, 41 appear twice, 3 appear 3 times
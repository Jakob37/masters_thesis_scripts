#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# Function input
# 
#
#

# (1) GET SIGNATURE VALUES  ---------------------------------------------------
# for proliferation or immune analysis

get_signature_value <- function(gene_symbol_column, sample_column, module, dup_removal_fun) {
  # returns a signature value for a matrix of one sample (no row names, two columns (gene symbols, sample), values are FPKM)
  # create a subset with gene symbols and sample
  current.sample <- bind_cols(data.frame(gene_symbol_column), data.frame(sample_column))
  # remove duplicate rows by getting max value for each symbol
  nodup <- current.sample %>% group_by(geneSymbol) %>% summarise_all(dup_removal_fun)
  # get only common genes
  only_common <- nodup %>% filter(nodup$geneSymbol %in% common_genes)
  # order by value
  ordered <- only_common %>% arrange(.[[2]])
  # extract position for genes in module, sum them and add to final table
  # fred
  signature.rows <- which(ordered$geneSymbol %in% module)
  signature.value <- sum(as.integer(signature.rows))
  return(signature.value)
}
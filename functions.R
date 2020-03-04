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


# (2) GET SAMPLES CLASS FROM AIMS RESULT  ---------------------------------------
# for CLAMS, ROR SSPs, PAM50 SSPs
get_samples_class_from_predictor_result <- function(aims_result) {
  sample.class <- aims_result$cl[,]
  sample.class <- data.frame(as.list(sample.class))
  sample.class <- t(sample.class)
  sample.class <- rownames_to_column(data.frame(sample.class), var = "sample.id")
}


# (3) GET SAMPLES PROB FROM AIMS RESULT  ---------------------------------------
# for CLAMS, ROR SSPs, PAM50 SSPs
get_samples_prob_from_predictor_result <- function(aims_result) {
  sample.prob <- aims_result$prob[,]
  sample.prob <- data.frame(as.list(sample.prob))
  sample.prob <- t(sample.prob)
  sample.prob <- rownames_to_column(data.frame(sample.prob), var = "sample.id")
}

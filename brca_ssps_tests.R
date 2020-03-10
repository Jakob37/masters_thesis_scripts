


# load things from data_input.R
# might need variables from dataset_report.R
# source functions.R



# General info on ROR SSPs  ----------------------------------

## Reduced set (ror.red.aims.gs)
# unique genes
red_genes <- get.unique.genes(ror.red.aims.gs[["all.pairs"]])
length(red_genes) # 492 unique genes
# rules
red_rules <- ror.red.aims.gs[["all.pairs"]]
length(red_rules) # 803 unique rules
# rules/gene
red_genes_in_rules <- get.rules.per.gene(ror.red.aims.gs[["all.pairs"]])
ggplot(mapping = aes(red_genes_in_rules)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1:19,29)) +
  labs(x = "Rules per gene", y ="Gene percentage") +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", 
            hjust=-0.2, vjust = 0.5, size=3, angle=90)
# 50 rules per class
red_classes <- names(ror.red.aims.gs[['selected.pairs.list']])
# maybe -----
# make a table of which genes are in which classes and how many times
# this would give an idea on how much missing genes will influence the analysis



## All 19k set (ror.all.aims.gs)
# unique genes
all_genes <- get.unique.genes(ror.all.aims.gs[["all.pairs"]])
length(all_genes) # 296 unique genes
# rules
all_rules <- ror.all.aims.gs[["all.pairs"]]
length(all_rules) # 361 unique rules
# rules/gene
all_genes_in_rules <- get.rules.per.gene(ror.all.aims.gs[["all.pairs"]])
ggplot(mapping = aes(all_genes_in_rules)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1:13,15,17)) +
  labs(x = "Rules per gene", y ="Gene percentage") +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", 
            angle=90, hjust=-0.2, vjust = 0.5, size=3)
# 21 rules per class
all_classes <- names(ror.all.aims.gs[['selected.pairs.list']])


# How many genes in common in the rules?
length(intersect(red_genes, all_genes)) # 163
# How many genes in rules in total, both models?
length(union(red_genes, all_genes)) # 625


### Which entrezIds are present in datasets

## Reduced set
# GOBO
sum(red_genes %in% entrez_ids_gobo) # 492/492, all

# SCAN-B
sum(red_genes %in% entrez_ids_scanb) # 489/492, 99.4%

# TCGA
sum(red_genes %in% entrez_ids_tcga) # 492/492, all



## All 19k set
# GOBO
sum(all_genes %in% entrez_ids_gobo) # 215/296, 72.6%

# SCAN-B
sum(all_genes %in% entrez_ids_scanb) # 286/296, 96.6%

# TCGA
sum(all_genes %in% entrez_ids_tcga) # 296/296, all





### Specific AIMS results


## Reduced set (ror_red_result_)

# GOBO
# no warnings raised

# SCAN-B
# missing/have more than one (3 rules): e10460<e9812, e10087<e29127, e3838<e6301

# TCGA
# no warnings raised


## All 19k set (ror_all_result_)

# GOBO
# missing/have more than one (153 rules): e55771<e79018,e23582<e6241,e10569<e9232,e113130<e9640,
#e29928<e9212,e150468<e9860,e890<e92345,e29089<e405,e29089<e63931,e11065<e84317,
#e153222<e699,e55165<e79699,e113130<e8924,e127253<e6790,e57658<e9232,e151648<e376940,
#e388630<e90381,e132884<e4288,e220929<e22974,e3680<e9232,e387119<e9787,e126669<e51203,
#e150468<e54413,e101060200<e55635,e116224<e890,e65062<e9833,e257364<e3833,e23013<e9232,
#e10403<e90826,e201292<e332,e10112<e57381,e353088<e3832,e24137<e6330,e100127983<e55355,
#e91653<e991,e6854<e90381,e4093<e55789,e56849<e64151,e55355<e6330,e1058<e256691,
#e25943<e699,e57685<e9833,e167691<e64151,e11065<e4629,e128553<e55771,e55165<e57381,
#e154810<e9133,e54756<e890,e5347<e91653,e151648<e8626,e84417<e9156,e25943<e55165,
#e126669<e55355,e6330<e8318,e1062<e285671,e4605<e83987,e113130<e57636,e127428<e150468,
#e151648<e56849,e55165<e57685,e167691<e9833,e11004<e6330,e389432<e55635,e146433<e4605,
#e256691<e699,e115825<e24137,e23005<e2305,e55165<e65062,e114928<e699,e154810<e22974,
#e6330<e9156,e389432<e9928,e10721<e80731,e146909<e388630,e4605<e91653,e57685<e9156,
#e64388<e90381,e332<e4629,e1058<e84417,e6330<e63967,e150468<e256691,e113115<e79633,
#e23397<e57381,e54510<e55165,e113130<e25802,e151648<e256691,e152503<e22974,e4605<e4629,
#e150468<e6330,e83540<e84986,e57685<e79801,e23235<e81610,e332<e64129,e55771<e57381,
#e55388<e6330,e81610<e9639,e113130<e257364,e151648<e79684,e26249<e90381,e54443<e80310,
#e373<e83540,e222236<e9156,e6330<e9133,e11065<e91653,e51203<e57381,e389072<e701,
#e4751<e57685,e202052<e24137,e100127983<e9833,e259266<e389432,e130497<e4605,e256691<e890,
#e116224<e51203,e22974<e6487,e100127983<e890,e202052<e699,e123879<e9133,e54443<e79817,
#e389072<e4288,e55247<e85480,e84449<e9493,e23397<e57685,e23594<e4093,e116159<e6241,
#e220134<e389432,e100127983<e9156,e10112<e6330,e5347<e56963,e4147<e9232,e27153<e4605,
#e10157<e150468,e11004<e51058,e151648<e5994,e116224<e9833,e24137<e64802,e51277<e55789,
#e1058<e115825,e119032<e983,e1058<e285268,e29089<e51232,e150468<e376940,e151648<e54756,
#e4040<e81610,e113130<e29907,e55696<e9232,e113115<e167691,e10580<e54443,e100506144<e9133,
#e113130<e23371,e150468<e26137,e113115<e1524,e1869<e84892,e1185<e151648
length(ror_all_result_gobo[["EntrezID.used"]])


# SCAN-B
# missing/have more than one (20 rules): e23371<e4605,e151648<e376940,e10460<e23371,
# e151648<e8626,e84417<e9156,e151648<e56849,e10186<e11065,e1058<e84417,e151648<e256691,
# e151648<e79684,e10087<e29127,e151648<e5994,e119032<e983,e23506<e890,e151648<e54756,
# e100506144<e9133,e220972<e23397,e113130<e23371,e332<e55793,e1185<e151648

# TCGA
# no warnings raised


# General info on PAM50 SSPs  ----------------------------------

## Reduced set (pam50.red.aims.gs)
# unique genes
red_genes <- get.unique.genes(pam50.red.aims.gs[["all.pairs"]])
length(red_genes) # 161 unique genes
# rules
red_rules <- pam50.red.aims.gs[["all.pairs"]]
length(red_rules) # 88 unique rules
# rules/gene
red_genes_in_rules <- get.rules.per.gene(pam50.red.aims.gs[["all.pairs"]])
ggplot(mapping = aes(red_genes_in_rules)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1:19,29)) +
  labs(x = "Rules per gene", y ="Gene percentage") +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", 
            hjust=-0.2, vjust = 0.5, size=3, angle=90)
# 91.3% of genes appear in only one rule, 8.1% in 2, 0.6% in 3.
# 22 rules per class (4 classes)
red_classes <- names(pam50.red.aims.gs[['selected.pairs.list']])
# maybe -----
# make a table of which genes are in which classes and how many times
# this would give an idea on how much missing genes will influence the analysis



## All 19k set (pam50.all.aims.gs)
# unique genes
all_genes <- get.unique.genes(pam50.all.aims.gs[["all.pairs"]])
length(all_genes) # 153 unique genes
# rules
all_rules <- pam50.all.aims.gs[["all.pairs"]]
length(all_rules) # 84 unique rules
# rules/gene
all_genes_in_rules <- get.rules.per.gene(pam50.all.aims.gs[["all.pairs"]])
ggplot(mapping = aes(all_genes_in_rules)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1:13,15,17)) +
  labs(x = "Rules per gene", y ="Gene percentage") +
  scale_y_continuous(labels=scales::percent) +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", 
            angle=90, hjust=-0.2, vjust = 0.5, size=3)
# 90% of genes in only 1 rule, 10% in 2 rules.
# 21 rules per class (4 classes)
all_classes <- names(pam50.all.aims.gs[['selected.pairs.list']])


# How many genes in common in the rules?
length(intersect(red_genes, all_genes)) # 94
# How many genes in rules in total, both models?
length(union(red_genes, all_genes)) # 220


### Which entrezIds are present in datasets

## Reduced set
# GOBO
sum(red_genes %in% entrez_ids_gobo) # 161/161, all

# SCAN-B
sum(red_genes %in% entrez_ids_scanb) # 161/161, all

# TCGA
sum(red_genes %in% entrez_ids_tcga) # 161/161, all



## All 19k set
# GOBO
sum(all_genes %in% entrez_ids_gobo) # 114/153, 74.51%

# SCAN-B
sum(all_genes %in% entrez_ids_scanb) # 150/153, 98.04%

# TCGA
sum(all_genes %in% entrez_ids_tcga) # 153/153, all





### Specific AIMS results


## Reduced set (pam50_red_result_)

# GOBO
# no warnings raised

# SCAN-B
# no warnings raised

# TCGA
# no warnings raised


## All 19k set (pam50_all_result_)

# GOBO
# missing/have more than one (39 rules): e51755<e85004,e10948<e55861,e50937<e84879,e23158<e84299,e144110<e5046,
#e2264<e400793,e285386<e79801,e167691<e64151,e388630<e90381,e116224<e890,e4605<e91653,e11065<e57727,e55355<e6330,
#e10403<e90826,e22974<e84181,e23331<e81610,e57381<e9833,e10283<e9232,e257364<e332,e128239<e23650,e147841<e3868,
#e151648<e2296,e342667<e79733,e5653<e90381,e10721<e81706,e150468<e5268,e1410<e9232,e29089<e3861,e222171<e2296,
#e319089<e53335,e1824<e5205,e155465<e2568,e79875<e84441,e25803<e81628,e83879<e8416,e375035<e79818,e29968<e401546,
#e388743<e81706,e140578<e1767

# SCAN-B
# missing/have more than one (3 rules): e333926<e79919,e151648<e2296,e319089<e53335

# TCGA
# no warnings raised


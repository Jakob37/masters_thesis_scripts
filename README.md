# Masters in Bioinformatics - thesis scripts

Hello Deborah!

Scripts created during my Masters thesis at Lund University: "Pan-cancer validation of a lung-adenocarcinoma-derived gene-expression-based prognostic predictor". It basically involves running the predictor named CLAMS on tumor samples from 32 different cancer types and analysing the results.  
Different analyses use different scripts, all written in RStudio running R version 3.6.3. Irrespect of what you want to analyze, you start by running parts of data_input.R to load the data.  
The only script intended to be sourced is functions.R. For all others, run code line by line or in chunks according to computer capacity.  
The following packages are needed in some or all scripts: Biobase v2.46.0, CLAMS (still not released to the public - Liljedahl et al., submitted), e1071 v1.7-3, org.Hs.eg.db v3.10.0, pheatmap v1.0.12, survminer v0.4.6, survival v3.1-11, tidyverse v1.3.0, viridisLite v0.3.0, yaml v2.2.1.  

## Configuration files

* config_alltcga.yml  
* config_brca_gobo.yml  
* config_brca_scanb.yml  

Configuration files for the three datasets (TCGA, GOBO, SCAN-B) with file paths, specific column names in tables, etc. Used mainly in previous versions of the scripts.  
Packages: none, but yaml to read them.

## Basic scripts

* data_input.R

Contains file paths to all files needed in other scripts, such as patient annotation tables, gene expression matrices, single sampre predictors (SSPs), classification tables after SSP analyis, etc.  
Packages: tidyverse, yaml, Biobase, org.Hs.eg.db.  

* functions.R

Functions needed in different scripts for getting gene rank value, classification from SSP results, count in how many rules a gene appears, etc.  
Packages: none.  

* plots_names.R

Information needed for changing axis information in plots. Contains cancer type acronyms and corresponding full names, acronyms with asterisks marking statistical significancy, etc.  
Packages: none.  

* run_ssps_function.R

Runs different SSPs, based on github.com/meoyo/trainAIMS. Used for ROR and PAM50 SSPs, CLAMS does not need it.  
Packages: e1071.  

## Analysis scripts

* clams_pancancer_analysis.R

Runs CLAMS on all samples from the three datasets. Check for congruency with other prognostic subtypes.  
Packages: tidyverse, yaml, org.Hs.eg.db, CLAMS, e1071.  

* OS.R

Overall survival analysis. Takes as input one table with information for all three datasets, mutate columns as needed (e.g., all to months or years), censor at any time point. Used with the following analyses: CLAMS, proliferation, ROR, PAM50.  
Packages: tidyverse, survival, survminer, viridisLite.  

* cox.R

Cox regression analysis. Just as overall survival, takes as input one table with information for all three datasets, mutate columns as needed (e.g., all to months or years), censor at any time point. Used with the following analyses: CLAMS, proliferation, ROR, PAM50.  
Packages: tidyverse, survival, survminer.  


* brca_pam50_pancancer.R

Runs the PAM50 SSP on all three datasets, compares information for two versions of this predictor (trained on a full gene set and on a reduced gene set).  
Packages: tidyverse, org.Hs.eg.db.  

* brca_ror_pancancer.R

Runs the ROR SSP on all three datasets, compares information for two versions of this predictor (trained on a full gene set and on a reduced gene set).  
Packages: tidyverse, org.Hs.eg.db, viridisLite.  

* proliferation_analysis.R

Uses genes associated with cell proliferation to show how proliferative tumor samples are. Compares these between CLAMS groups. Values used to split samples into a low- and a high-proliferative group.  
Packages: tidyverse.  

* immune_analysis.R

Similar to proliferation_analysis.R, uses an immune signature to see if tumor samples vary.  
Packages: tidyverse.  


## Other

* brca_ssps_tests.R

Gets general information from the ROR and PAM50 SSPs and the analysis result for each, such as how many genes are used in the rules and how many rules were missing in the analysis.  
Packages: none.  

* dataset_report.R

Gets basic information from TCGA, GOBO and SCAN-B, such as how many duplicated rows there are in the gene expression matrix.  
Packages: tidyverse, org.Hs.eg.db, viridisLite.  

* ssps_heatmap.R

Creates heatmaps (samples x rules) for CLAMS, PAM50, and ROR.  
Packages: tidyverse, pheatmap, viridisLite, CLAMS.

* summary_analysis.R

Creates cross-dataset summary plots of sample sizes and CLAMS distribution.  
Packages: tidyverse.   

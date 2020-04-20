# Masters in Bioinformatics - thesis scripts

Scripts created during my Masters thesis to run CLAMS (prognosis predictor developed for lung adenocarcinoma) pan-cancer in R.  
Different analyses use different scripts. Irrespect of what you want to analyze, you start by running parts of data_input.R to load the data.  
The only script intended to be sourced is functions.R. For all others, run code line by line or in chunks according to computer capacity.  

## Configuration files

* config_alltcga.yml  
* config_brca_gobo.yml  
* config_brca_scanb.yml  

Configuration files for the three datasets (TCGA, GOBO, SCAN-B) with file paths, specific column names in tables, etc. Used mainly in previous versions of the scripts.

## Basic scripts

* data_input.R

Contains file paths to all files needed in other scripts, such as patient annotation tables, gene expression matrices, single sampre predictors (SSPs), classification tables after SSP analyis, etc.

* functions.R

Functions needed in different scripts for getting gene rank value, classification from SSP results, count in how many rules a gene appears, etc.

* plots_names.R

Information needed for changing axis information in plots. Contains cancer type acronyms and corresponding full names, acronyms with asterisks marking statistical significancy, etc.

* run_ssps_function.R

Runs different SSPs, based on github.com/meoyo/trainAIMS. Used for ROR and PAM50 SSPs, CLAMS does not need it.

## Analysis scripts

* clams_pancancer_analysis.R

Runs CLAMS on all samples from the three datasets. Check for congruency with other prognostic subtypes.

* OS.R

Overall survival analysis. Takes as input one table with information for all three datasets, mutate columns as needed (e.g., all to months or years), censor at any time point. Used with the following analyses: CLAMS, proliferation, ROR, PAM50.

* cox.R

Cox regression analysis. Just as overall survival, takes as input one table with information for all three datasets, mutate columns as needed (e.g., all to months or years), censor at any time point. Used with the following analyses: CLAMS, proliferation, ROR, PAM50.


* brca_pam50_pancancer.R

Runs the PAM50 SSP on all three datasets, compares information for two versions of this predictor (trained on a full gene set and on a reduced gene set).

* brca_ror_pancancer.R

Runs the ROR SSP on all three datasets, compares information for two versions of this predictor (trained on a full gene set and on a reduced gene set).


* proliferation_analysis.R

Uses genes associated with cell proliferation to show how proliferative tumor samples are. Compares these between CLAMS groups. Values used to split samples into a low- and a high-proliferative group.

* immune_analysis.R

Similar to proliferation_analysis.R, uses an immune signature to see if tumor samples vary.


## Other

* brca_ssps_tests.R

Gets general information from the ROR and PAM50 SSPs and the analysis result for each, such as how many genes are used in the rules and how many rules were missing in the analysis.  

* dataset_report.R

Gets basic information from TCGA, GOBO and SCAN-B, such as how many duplicated rows there are in the gene expression matrix.

* ssps_heatmap.R

Creates heatmaps (samples x rules) for CLAMS, PAM50, and ROR.

* summary_analysis.R

Creates cross-dataset summary plots of sample sizes and CLAMS distribution.

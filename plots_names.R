#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS
#
# Cancer names for plot legends
#
#



full_names <- c('THCA_TCGA' = 'Thyroid carcinoma (THCA)', 
                'KIRP_TCGA' = 'Kidney renal papillary cell carcinoma (KIRP)', 
                'PRAD_TCGA' = 'Prostate adenocarcinoma (PRAD)',
                'LUAD_TCGA' = 'Lung adenocarcinoma (LUAD)', 
                'KICH_TCGA' = 'Kidney chromophobe (KICH)', 
                'KIRC_TCGA' = 'Kidney renal clear cell carcinoma (KIRC)',
                'PCPG_TCGA' = 'Pheochromocytoma and paraganglioma (PCPG)', 
                'LGG_TCGA' = 'Brain lower grade glioma (LGG)', 
                'ACC_TCGA' = 'Adrenocortical carcinoma (ACC)',
                'LIHC_TCGA' = 'Liver hepatocellular carcinoma (LIHC)', 
                'PAAD_TCGA' = 'Pancreatic adenocarcinoma (PAAD)', 
                'BRCA_TCGA' = 'Breast invasive carcinoma (BRCA)',
                'MESO_TCGA' = 'Mesothelioma (MESO)', 
                'CHOL_TCGA' = 'Cholangiocarcinoma (CHOL)', 
                'THYM_TCGA' = 'Thymoma (THYM)', 
                'STAD_TCGA' = 'Stomach adenocarcinoma (STAD)',
                'LUSC_TCGA' = 'Lung squamous cell carcinoma (LUSC)', 
                'UVM_TCGA' = 'Uveal melanoma (UVM)', 
                'SARC_TCGA' = 'Sarcoma (SARC)', 
                'BLCA_TCGA' = 'Bladder urothelial carcinoma (BLCA)',
                'OV_TCGA' = 'Ovarian serous cystadenocarcinoma (OV)', 
                'GBM_TCGA' = 'Glioblastoma multiforme (GBM)', 
                'HNSC_TCGA' = 'Head and neck squamous cell carcinoma (HNSC)',
                'ESCA_TCGA' = 'Esophageal carcinoma (ESCA)', 
                'COAD_TCGA' = 'Colon adenocarcinoma (COAD)', 
                'READ_TCGA' = 'Rectum adenocarcinoma (READ)',
                'UCEC_TCGA' = 'Uterine corpus endometrial carcinoma (UCEC)', 
                'UCS_TCGA' = 'Uterine carcinosarcoma (UCS)',
                'CESC_TCGA' = 'Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)',
                'TGCT_TCGA' = 'Testicular germ cell tumors (TGCT)', 
                'SKCM_TCGA' = 'Skin cutaneous melanoma (SKCM)', 
                'DLBC_TCGA' = 'Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)',
                'BRCA_GOBO' = '[GOBO] Breast invasive carcinoma (BRCA)', 
                'BRCA_SCAN-B' = '[SCAN-B] Breast invasive carcinoma (BRCA)')

# graph_name_no_ast
acronym_only <- c('THCA_TCGA' = 'THCA', 'KIRP_TCGA' = 'KIRP', 'PRAD_TCGA' = 'PRAD',
                  'LUAD_TCGA' = 'LUAD', 'KICH_TCGA' = 'KICH', 'KIRC_TCGA' = 'KIRC',
                  'PCPG_TCGA' = 'PCPG', 'LGG_TCGA' = 'LGG', 'ACC_TCGA' = 'ACC',
                  'LIHC_TCGA' = 'LIHC', 'PAAD_TCGA' = 'PAAD', 'BRCA_TCGA' = 'BRCA',
                  'MESO_TCGA' = 'MESO', 'CHOL_TCGA' = 'CHOL', 'THYM_TCGA' = 'THYM', 'STAD_TCGA' = 'STAD',
                  'LUSC_TCGA' = 'LUSC', 'UVM_TCGA' = 'UVM', 'SARC_TCGA' = 'SARC', 'BLCA_TCGA' = 'BLCA',
                  'OV_TCGA' = 'OV', 'GBM_TCGA' = 'GBM', 'HNSC_TCGA' = 'HNSC',
                  'ESCA_TCGA' = 'ESCA', 'COAD_TCGA' = 'COAD', 'READ_TCGA' = 'READ',
                  'UCEC_TCGA' = 'UCEC', 'UCS_TCGA' = 'UCS', 'CESC_TCGA' = 'CESC',
                  'TGCT_TCGA' = 'TGCT', 'SKCM_TCGA' = 'SKCM', 'DLBC_TCGA' = 'DLBC',
                  'BRCA_GOBO' = '[GOBO] BRCA', 'BRCA_SCAN-B' = '[SCAN-B] BRCA')

# correct_graph_name
acronym_ast_os_clams <- c('THCA_TCGA' = 'THCA', 'KIRP_TCGA' = '*KIRP', 'PRAD_TCGA' = 'PRAD',
                        'LUAD_TCGA' = '*LUAD', 'KICH_TCGA' = 'KICH', 'KIRC_TCGA' = '*KIRC',
                        'PCPG_TCGA' = 'PCPG', 'LGG_TCGA' = '*LGG', 'ACC_TCGA' = 'ACC',
                        'LIHC_TCGA' = '*LIHC', 'PAAD_TCGA' = 'PAAD', 'BRCA_TCGA' = 'BRCA',
                        'MESO_TCGA' = 'MESO', 'CHOL_TCGA' = 'CHOL', 'THYM_TCGA' = 'THYM', 'STAD_TCGA' = 'STAD',
                        'LUSC_TCGA' = 'LUSC', 'UVM_TCGA' = 'UVM', 'SARC_TCGA' = 'SARC', 'BLCA_TCGA' = 'BLCA',
                        'OV_TCGA' = 'OV', 'GBM_TCGA' = 'GBM', 'HNSC_TCGA' = 'HNSC',
                        'ESCA_TCGA' = 'ESCA', 'COAD_TCGA' = 'COAD', 'READ_TCGA' = 'READ',
                        'UCEC_TCGA' = 'UCEC', 'UCS_TCGA' = 'UCS', 'CESC_TCGA' = 'CESC',
                        'TGCT_TCGA' = 'TGCT', 'SKCM_TCGA' = 'SKCM', 'DLBC_TCGA' = 'DLBC',
                        'BRCA_GOBO' = '[GOBO] *BRCA', 'BRCA_SCAN-B' = '[SCAN-B] *BRCA')

#sign_os_prolif
acronym_ast_os_prolif_highlow <- c('THCA_TCGA' = 'THCA', 'KIRP_TCGA' = '*KIRP', 'PRAD_TCGA' = 'PRAD',
                                  'LUAD_TCGA' = '*LUAD', 'KICH_TCGA' = 'KICH', 'KIRC_TCGA' = 'KIRC',
                                  'PCPG_TCGA' = 'PCPG', 'LGG_TCGA' = '*LGG', 'ACC_TCGA' = '*ACC',
                                  'LIHC_TCGA' = '*LIHC', 'PAAD_TCGA' = '*PAAD', 'BRCA_TCGA' = 'BRCA',
                                  'MESO_TCGA' = '*MESO', 'CHOL_TCGA' = '*CHOL', 'THYM_TCGA' = 'THYM', 'STAD_TCGA' = 'STAD',
                                  'LUSC_TCGA' = 'LUSC', 'UVM_TCGA' = 'UVM', 'SARC_TCGA' = '*SARC', 'BLCA_TCGA' = 'BLCA',
                                  'OV_TCGA' = 'OV', 'GBM_TCGA' = 'GBM', 'HNSC_TCGA' = 'HNSC',
                                  'ESCA_TCGA' = 'ESCA', 'COAD_TCGA' = 'COAD', 'READ_TCGA' = 'READ',
                                  'UCEC_TCGA' = 'UCEC', 'UCS_TCGA' = 'UCS', 'CESC_TCGA' = 'CESC',
                                  'TGCT_TCGA' = 'TGCT', 'SKCM_TCGA' = 'SKCM', 'DLBC_TCGA' = 'DLBC',
                                  'BRCA_GOBO' = '[GOBO] *BRCA', 'BRCA_SCAN-B' = '[SCAN-B] *BRCA')
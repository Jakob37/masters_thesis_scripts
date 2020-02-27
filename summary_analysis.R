############################
#
# DEBORAH FIGUEIREDO NACER DE OLIVEIRA
# 2019-2020
#
# MASTER OF BIOINFORMATICS - THESIS 
#
# SCRIPT TO GET SUMMARY DATA AND CREATE CROSS-DATA PLOTS
#
############################

# CLEAR ENVIRONMENT ---------------------------------------------------------------------------------------------------

rm(list=ls())
while (dev.cur()>1) dev.off()


# LOAD GENERAL PACKAGES ---------------------------------------------------------------------------------------------------

library(tidyverse)


# INPUT DATA FROM data_input.R  ----------------------------------

file.edit("/media/deboraholi/Data/LUND/9 THESIS/src/data_input.R")

# CLAMS RESULTS 
# PLOTS NAMES


# PLOTS  ---------------------------------------------------------

patient_annotation_clams <- patient_annotation_clams %>% mutate(groups.to.analyze = paste(cancer.type, dataset, sep = '_'))

## Sample sizes

# barplot
patient_annotation_clams %>% 
  group_by(groups.to.analyze) %>% 
  tally() %>%
  ggplot(aes(x=reorder(groups.to.analyze, -n), y=n)) +
  geom_bar(stat="identity", fill="darkgrey") +
  coord_flip() +
  labs(x ="", y = "Number of samples") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.y = element_line(color="grey50"),
        axis.line.x = element_line(color="grey50"), axis.ticks.x.bottom=element_line(color="grey50"),
        legend.position="none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = full_names)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/0_clams/samples_per_type.png", width=9, height=6, dpi=300)


# lollipop
patient_annotation_clams %>% 
  group_by(groups.to.analyze) %>% 
  tally() %>%
  ggplot(aes(x=reorder(groups.to.analyze, -n), y=n)) +
  geom_segment( aes(x=reorder(groups.to.analyze, -n), 
                    xend=reorder(groups.to.analyze, -n), 
                    y=0, yend=n), color="darkgrey") +
  geom_point(size=2, color="darkgrey") +
  coord_flip() +
  labs(x ="", y = "Number of samples") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.y = element_line(color="grey50"),
        axis.line.x = element_line(color="grey50"), axis.ticks.x.bottom=element_line(color="grey50"),
        legend.position="none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3600)) +
  scale_x_discrete(labels = full_names)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/0_clams/samples_per_type_lollipop.png", width=9, height=6, dpi=300)


## Sample sizes and CLAMS

patient_annotation_clams %>% 
  group_by(groups.to.analyze, clams.class) %>% 
  tally() %>%
  pivot_wider(id_cols = groups.to.analyze, names_from = clams.class, values_from = n) %>% 
  mutate_at("TRU", ~replace(., is.na(.), 0))  %>% 
  mutate(tru_perc = (TRU / (TRU+NonTRU))) %>% mutate(total_samples = (TRU+NonTRU)) %>%
  arrange(desc(tru_perc)) %>% pivot_longer(cols=c("NonTRU", "TRU"), names_to = "clams.class") %>% 
  ggplot(aes(x=reorder(groups.to.analyze, tru_perc), y=value)) +
  geom_bar(position="fill", stat="identity", aes(fill=clams.class)) +
  theme_minimal() +
  coord_flip(clip = "off") +
  labs(x ="", y = "Sample percentage") +
  geom_text(aes(x=groups.to.analyze, y=1.06, label=total_samples), size=3, color="grey50") +
  geom_text(aes(x=34, y=1.13, label="= n"), size=3, color="grey50") +
  theme(axis.text.x = element_text(size = 11),
        panel.grid.major.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x.bottom=element_line(color="grey50")) +
  scale_fill_manual(values=c("deepskyblue3", "darkorange1"), guide=guide_legend(reverse = TRUE), name="CLAMS") +
  scale_x_discrete(labels = full_names) +
  scale_y_continuous(expand = c(0, 0), labels=scales::percent)

ggsave("/media/deboraholi/Data/LUND/9 THESIS/0_clams/clams_per_type.png", width=9, height=6, dpi=300)

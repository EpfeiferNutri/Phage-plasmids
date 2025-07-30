# libraries

library(data.table)
library(seqinr)
library(scales)
library(ggbreak)
library(writexl)
library(patchwork)
library("rhmmer")
library(readxl)
library(tidyverse)


# main table
All_MGE_information_table = read_excel("Study_associated_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

################# Curation of cp32-related P-Ps into a well-defined groups

wGRR_df = read_tsv("Study_associated_data/Supplementary_data/tyPPing_criteria_tables/cp32_0312_wGRR_df.tsv") %>% 
  mutate(replicon=query)

cp32_pps_0321 = wGRR_df %>% filter(subject %in% c("NC_012158.1", "NC_012156.1", "NC_012154.1", "NC_012152.1", "NC_012149.1"),
                                   !replicon %in% c("NC_012158.1", "NC_012156.1", "NC_012154.1", "NC_012152.1", "NC_012149.1")) %>% 
  filter(wGRR >= 0.0, replicon != subject, n_similar_proteins >= prot_subject/2 ) %>% 
  # count the hits per subject replicon
  group_by(replicon) %>% arrange(desc(wGRR)) %>% mutate(n_hits = 1:n()) %>%
  #just take the best hit (highest wGRR)
  filter(n_hits == 1)

################# Figure (A)
ggplot(cp32_pps_0321, aes(x = wGRR)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "wGRR value distribution of 03/21 plasmids from Borrelia and Borreliella species", x = "wGRR value", y = "Number of genomes") +
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme_minimal()+
  theme(
    strip.background = element_blank(),
    text = element_text(size = 20))

ref_cp32_pps_0321 = cp32_pps_0321 %>% filter(wGRR >= 0.6) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `Replicon size (bp)`), by=c("replicon"="NCBI accession"))

################# Figure (B)
ref_cp32_pps_0321 %>%
  ggplot(aes(x = `Replicon size (bp)` / 1000)) +
  geom_density(alpha = 0.8, fill='darkgreen', col = 'darkgreen') +
  theme_bw() +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  scale_x_continuous(breaks = seq(0, 60, 5)) +
  xlab("Genome size (kb)") +
  ylab("Density") +
  theme(
    strip.background = element_blank(),
    text = element_text(size = 22))

################################################



################# P-P type specific size distribution and ranges

# size distribution figure

PP_size = All_MGE_information_table %>%
  filter(`in 03/21` == "Yes", `MGE type MM-GRC` %in% c("AB_1 subgroup", "P1_1 subgroup", "P1_2 subgroup", "N15 group", "SSU5_pHCM2 group", 
                                         "pMT1 group", "pCAV group", "pSLy3 group", "pKpn group", "cp32 group"))

PP_known_size = PP_size %>% select(`NCBI accession`, `Replicon size (bp)`, `MGE type MM-GRC`) %>%
  mutate(`MGE type MM-GRC` = str_remove(`MGE type MM-GRC`, " (group|subgroup)$"))

peak_x_values = PP_known_size %>% group_by(`MGE type MM-GRC`) %>%
  summarise(peak_x = {
    dens = density(`Replicon size (bp)`)              # Calculate density
    dens$x[which.max(dens$y)]                  # Extract x at max density
  }, .groups = "drop") %>% mutate(peak_x = round(peak_x / 1000))

facet_counts = PP_known_size %>%
  group_by(`MGE type MM-GRC`) %>%
  summarise(count = n(), .groups = "drop")

# Merge peak_x_values with facet_counts for annotations
annotations = facet_counts %>%
  left_join(peak_x_values, by = "MGE type MM-GRC")

# Genome size distributions per P-P type 
PP_known_size %>%
  ggplot(aes(x = `Replicon size (bp)` / 1000)) +
  geom_density(alpha = 0.9, color = "darkgreen", fill='limegreen') +
  theme_bw() +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = seq(0, 230, 40)) +
  facet_wrap(~`MGE type MM-GRC`, scales = "free_y", nrow = 2) +
  xlab("Genome size by P-P type (kb)") +
  ylab("Density") +
  geom_text(
    data = annotations,
    aes(x = Inf, y = Inf, label = paste("N =", count, "\nPeak = ", round(peak_x, 1), " kb")),
    inherit.aes = FALSE,
    hjust = 1.1,  
    vjust = 1.6, 
    size = 5, color = "black"
  ) +
  theme(
    strip.background = element_blank(),
    text = element_text(size = 22)
  )


################# Cutoffs for MinProteins 

MinGenes_cutoff = read_excel("Study_associated_data/Supplementary_data/S2_Profile_information_table.xlsx")  %>% 
  select(hmm_profile = `HMM profile ID`, sequence_score_thr = `Sequence score threshold`)

### proten IDs from 0321 are kept
protein_to_genome = fread("Study_associated_data/Supplementary_data/tyPPing_criteria_tables/protein_to_genome_0321.tsv")

### proten IDs from 0321 are used in HMM_search_output; 
HMM_search_output = read_domtblout('Study_associated_data/Supplementary_data/tyPPing_criteria_tables/all_pers_hmm_phages_plasmids_0321.tbl.out') %>%
  mutate(ali_length = abs(ali_to - ali_from+1), cov_profile = round(ali_length/qlen,3)) %>% 
  select(protein_id = domain_name, hmm_profile = query_name, domain_ievalue, cov_profile, sequence_score) %>%
  inner_join(protein_to_genome %>% select(protein_id = `protein_id 0321`, genome_id), by="protein_id") %>%
  # count the hits per genome_id
  group_by(genome_id, hmm_profile) %>% 
  arrange(domain_ievalue) %>% 
  mutate(n_hits_profile = 1:n()) %>%
  #just take the best hit (lowest domain_ievalue)
  filter(n_hits_profile == 1) %>%
  mutate(PP_category = str_extract(hmm_profile, pattern=".*(?=_pers_)"))

# HMM_search_output$domain_name %>% head()

HMM_filtered_MinGenes = left_join(HMM_search_output, MinGenes_cutoff, by="hmm_profile") %>% filter(sequence_score >= sequence_score_thr)

### AB_1

HMM_filtered_AB = HMM_filtered_MinGenes %>% filter(PP_category == "AB")  

phages_plasmids_0321_AB_hmm_processed = HMM_filtered_AB %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "AB_1 subgroup" ~ "AB_1 subgroup", TRUE ~ "P-P (other)"))

phages_plasmids_0321_AB_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 50, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 100, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("AB_1 subgroup" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### N15

HMM_filtered_N15 = HMM_filtered_MinGenes %>% filter(PP_category == "N15")  

phages_plasmids_0321_N15_hmm_processed = HMM_filtered_N15 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "N15 group" ~ "N15 group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_N15_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 40, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 100, 2)) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_fill_manual(values = c("N15 group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### P1_1

HMM_filtered_P1_1 = HMM_filtered_MinGenes %>% filter(PP_category == "P11")  

phages_plasmids_0321_P1_1_hmm_processed = HMM_filtered_P1_1 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "P1_1 subgroup" ~ "P1_1 subgroup", TRUE ~ "P-P (other)"))

phages_plasmids_0321_P1_1_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 40, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 130, 10)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("P1_1 subgroup" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  scale_y_break(c(44, 100)) +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### P1_2

HMM_filtered_P1_2 = HMM_filtered_MinGenes %>% filter(PP_category == "P12")  

phages_plasmids_0321_P1_2_hmm_processed = HMM_filtered_P1_2 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "P1_2 subgroup" ~ "P1_2 subgroup", TRUE ~ "P-P (other)"))

phages_plasmids_0321_P1_2_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 40, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("P1_2 subgroup" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  scale_y_break(c(15, 30)) + scale_y_break(c(40,185)) +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### cp32

HMM_filtered_cp32 = HMM_filtered_MinGenes %>% filter(PP_category == "cp32")  

phages_plasmids_0321_cp32_hmm_processed = HMM_filtered_cp32 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "cp32 group" ~ "cp32 group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_cp32_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 28, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 160, 5)) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_fill_manual(values = c("cp32 group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### pCAV

HMM_filtered_pCAV = HMM_filtered_MinGenes %>% filter(PP_category == "pCAV")  

phages_plasmids_0321_pCAV_hmm_processed = HMM_filtered_pCAV %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "pCAV group" ~ "pCAV group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_pCAV_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 40, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 160, 1)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("pCAV group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### pMT1

HMM_filtered_pMT1 = HMM_filtered_MinGenes %>% filter(PP_category == "pMT1")  

phages_plasmids_0321_pMT1_hmm_processed = HMM_filtered_pMT1 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "pMT1 group" ~ "pMT1 group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_pMT1_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 40, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 700, 10)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("pMT1 group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +  scale_y_break(c(35, 280)) +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### pSLy3

HMM_filtered_pSLy3 = HMM_filtered_MinGenes %>% filter(PP_category == "pSLy3")  

phages_plasmids_0321_pSLy3_hmm_processed = HMM_filtered_pSLy3 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "pSLy3 group" ~ "pSLy3 group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_pSLy3_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 60, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 700, 4)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("pSLy3 group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +  scale_y_break(c(30, 60)) +
  theme(strip.background = element_blank(),text = element_text(size = 25))


### SSU5_pHCM2

HMM_filtered_SSU5_pHCM2 = HMM_filtered_MinGenes %>% filter(PP_category == "SSU5")  

phages_plasmids_0321_SSU5_pHCM2_hmm_processed = HMM_filtered_SSU5_pHCM2 %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "SSU5_pHCM2 group" ~ "SSU5_pHCM2 group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_SSU5_pHCM2_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 60, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 700, 5)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("SSU5_pHCM2 group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +  
  theme(strip.background = element_blank(),text = element_text(size = 25))



### pKpn

HMM_filtered_pKpn = HMM_filtered_MinGenes %>% filter(PP_category == "pKpn")  

phages_plasmids_0321_pKpn_hmm_processed = HMM_filtered_pKpn %>% 
  group_by(genome_id) %>% summarise(n_profiles_per_replicon = n()) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`), by=c("genome_id"="NCBI accession")) %>%
  mutate(`MGE type` = case_when(`MGE type MM-GRC` == "Plasmid" ~ "Plasmid", `MGE type MM-GRC` == "Phage" ~ "Phage", 
                                `MGE type MM-GRC` == "pKpn group" ~ "pKpn group", TRUE ~ "P-P (other)"))

phages_plasmids_0321_pKpn_hmm_processed %>%
  ggplot(aes(x = n_profiles_per_replicon, fill = `MGE type`)) +
  geom_histogram(bins = 60, col = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 700, 4)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("pKpn group" = "green", "Plasmid" = "orange", "Phage" = "cornflowerblue", "P-P (other)" = 'gray')) +
  xlab("# Profile matches per genome") +
  ylab("Counts") +  
  theme(strip.background = element_blank(),text = element_text(size = 25))




################# Search for cutoffs for profile-to-protein hits (example P1-like sequences) for Composition

# P1_1

Table_P1_1_test = fread("Study_associated_data/Supplementary_data/tyPPing_criteria_tables/Table_P1_1_composition_positive_hit_cutoff.tsv")
Table_P1_1_test$`Positive hit cutoff` <- factor(
  Table_P1_1_test$`Positive hit cutoff`,
  levels = c(
    "No cutoffs",
    "Coverage >50%",
    "Coverage >50%, E-value < 1e-3",
    "Coverage >50%, E-value < 1e-20"))

# overview plot
Table_P1_1_test %>% ggplot(aes(x=`Tolerance to gaps`, y = `Number of genomes`, fill = `MGE type`)) + 
  geom_col(width = 0.8)+
  scale_fill_manual(values = c("P1_1" = "limegreen", "Plasmid" = "orange", "P1_2" = "gray", "P-P (other)" = 'gray50')) +
  facet_wrap(~`Positive hit cutoff`) +
  scale_x_continuous(breaks = seq(0,25,3), limits = c(-1,26) ,name = "# Allowed missing signature proteins")+
  scale_y_continuous(breaks = seq(0,340,50), name = "# Detected genomes")+
  theme_bw()+theme(text = element_text(size = 28))+coord_flip()



################# Search for sizes of the P-P type specific composition including tolerances for miss-matches (example P1-like sequences)

# P1_1

Table_P1_1_test = fread("Study_associated_data/Supplementary_data/tyPPing_criteria_tables/Table_P1_1_composition_pattern_size_cutoff.tsv")
Table_P1_1_test$`Minimum % of pattern size` <- factor(
  Table_P1_1_test$`Minimum % of pattern size`,
  levels = c(
    "Any pattern size", 
    ">25% of pattern size",
    ">50% of pattern size",
    ">75% of pattern size",
    ">80% of pattern size",
    ">= than minimum pattern size"))

# overview plot
Table_P1_1_test %>% ggplot(aes(x=`Tolerance to gaps`, y = `Number of genomes`, fill = `MGE type`)) + 
  geom_col(width = 0.8)+
  scale_fill_manual(values = c("P1_1" = "limegreen", "Plasmid" = "orange", "P1_2" = "gray", "P-P (other)" = 'gray50')) +
  facet_wrap(~`Minimum % of pattern size`) +
  scale_x_continuous(breaks = seq(0,25,3), limits = c(-1,26) ,name = "# Allowed missing signature proteins")+
  scale_y_continuous(breaks = seq(0,340,50), name = "# Detected genomes")+
  theme_bw()+theme(text = element_text(size = 25))+coord_flip()



################# Chosen parameters for Composition

#### All types, 50% cov, 75% pattern

Table_all_composition_tolerance_to_gaps_cutoff = fread("Study_associated_data/Supplementary_data/tyPPing_criteria_tables/Table_all_composition_tolerance_to_gaps_cutoff.tsv") 

# overview plot
Table_all_composition_tolerance_to_gaps_cutoff %>%
  ggplot(aes(x = `Tolerance to gaps`, y = `Number of genomes`, fill = fill_color)) +
  geom_col(width = 0.8) +
  facet_wrap(~`P-P type tested`, scales = "free_x", nrow=2) + 
  scale_fill_identity(
    guide = "legend", 
    labels = c("P-P type (correct)", "Plasmid", "Phage", "P-P (other)"), 
    breaks = c("limegreen", "orange", "cornflowerblue", "gray50"),
    name = "MGE type"
  ) +
  scale_x_continuous(name = "# Allowed missing signature proteins", breaks = pretty_breaks(n = 10), limits = c(-1,19)) +
  scale_y_continuous(name = "# Detected genomes", breaks = pretty_breaks(n = 6)) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  coord_flip()


###############################################################################################




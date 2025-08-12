# libraries

library(data.table)
library(seqinr)
library(scales)
library(ggbreak)
library(patchwork)
library("rhmmer")
library(readxl)
library(tidyverse)


# main table
All_MGE_information_table = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

################# Figure SX1: Curation of cp32 community

wGRR_df = read_tsv("tyPPing_criteria_tables/cp32_0312_wGRR_df.tsv") %>% 
  mutate(replicon=query)

cp32_pps_0321 = wGRR_df %>% filter(subject %in% c("NC_012158.1", "NC_012156.1", "NC_012154.1", "NC_012152.1", "NC_012149.1"),
                                   !replicon %in% c("NC_012158.1", "NC_012156.1", "NC_012154.1", "NC_012152.1", "NC_012149.1")) %>% 
  filter(wGRR >= 0.0, replicon != subject, n_similar_proteins >= prot_subject/2 ) %>% 
  # count the hits per subject replicon
  group_by(replicon) %>% arrange(desc(wGRR)) %>% mutate(n_hits = 1:n()) %>%
  #just take the best hit (highest wGRR)
  filter(n_hits == 1)

################# Figure SX1 (A)
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

################# Figure SX1 (B)
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



################# Figure SX6: Genome size distributions by P-P type

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



################################################################################################

MinProteins_cutoff = read_excel("Paper_related_data/Supplementary_data/S2_Profile_information_table.xlsx")  %>% 
  select(hmm_profile = `HMM profile ID`, sequence_score_thr = `Sequence score threshold`)

### proten IDs from 0321 are kept
protein_to_genome = fread("tyPPing_criteria_tables/protein_to_genome_0321.tsv")

### proten IDs from 0321 are used in HMM_search_output; 
HMM_search_output = read_domtblout('tyPPing_criteria_tables/all_pers_hmm_phages_plasmids_0321.tbl.out') %>%
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

HMM_filtered_MinProteins = left_join(HMM_search_output, MinProteins_cutoff, by="hmm_profile") %>% filter(sequence_score >= sequence_score_thr)


################################################


### AB

HMM_filtered_AB = HMM_filtered_MinProteins %>% filter(PP_category == "AB")  

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

HMM_filtered_N15 = HMM_filtered_MinProteins %>% filter(PP_category == "N15")  

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

HMM_filtered_P1_1 = HMM_filtered_MinProteins %>% filter(PP_category == "P11")  

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

HMM_filtered_P1_2 = HMM_filtered_MinProteins %>% filter(PP_category == "P12")  

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

HMM_filtered_cp32 = HMM_filtered_MinProteins %>% filter(PP_category == "cp32")  

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

HMM_filtered_pCAV = HMM_filtered_MinProteins %>% filter(PP_category == "pCAV")  

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

HMM_filtered_pMT1 = HMM_filtered_MinProteins %>% filter(PP_category == "pMT1")  

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

HMM_filtered_pSLy3 = HMM_filtered_MinProteins %>% filter(PP_category == "pSLy3")  

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

HMM_filtered_SSU5_pHCM2 = HMM_filtered_MinProteins %>% filter(PP_category == "SSU5")  

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

HMM_filtered_pKpn = HMM_filtered_MinProteins %>% filter(PP_category == "pKpn")  

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




################################################################################################

##### Positive hit filter for Composition branch

MinProteins_cutoff = read_excel("Paper_related_data/Supplementary_data/S2_Profile_information_table.xlsx")  %>% 
  select(hmm_profile = `HMM profile ID`, sequence_score_thr = `Sequence score threshold`)

### proten IDs for 03/21 dataset
protein_to_genome = fread("tyPPing_criteria_tables/protein_to_genome_0321.tsv")

### HMM_search_output for 03/21 dataset; 
HMM_search_output = read_domtblout('tyPPing_criteria_tables/all_pers_hmm_phages_plasmids_0321.tbl.out') %>%
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

hits_MinProteins = left_join(HMM_search_output, MinProteins_cutoff, by="hmm_profile") %>% filter(sequence_score >= sequence_score_thr) %>% mutate(hit_thr = "MinProteins")
hits_cov50 = HMM_search_output %>% filter(cov_profile >= 0.5) %>% mutate(hit_thr = "Coverage >= 50%")
hits_cov50_evalue3 = HMM_search_output %>% filter(cov_profile >= 0.5, domain_ievalue < 1e-3) %>% mutate(hit_thr = "Coverage >= 50% & E-value < 1e-3")
hits_cov50_evalue20 = HMM_search_output %>% filter(cov_profile >= 0.5, domain_ievalue < 1e-20) %>% mutate(hit_thr = "Coverage >= 50% & E-value < 1e-20")
hits_any = HMM_search_output %>% mutate(hit_thr = "No cutoff")

hits_all = bind_rows(hits_any, hits_cov50, hits_cov50_evalue3, hits_cov50_evalue20, hits_MinProteins, ) %>% 
  left_join(All_MGE_information_table %>% select(`NCBI accession`, `MGE type MM-GRC`) %>% mutate(`MGE type MM-GRC` = case_when(
    `MGE type MM-GRC` == "Plasmid" ~ "Plasmid",
    `MGE type MM-GRC` == "Phage" ~ "Phage",
    TRUE ~ "P-P")), by=c("genome_id"="NCBI accession"))

# Reorder hit_thr by total hits
hits_all = hits_all %>% mutate(hit_thr = reorder(hit_thr, -table(hit_thr)[hit_thr]))


# Create a summary table with counts
summary_df = tibble(
  `Positive hit filter` = c(
    "Sequence score \n(as in MinProteins)\n",
    "Coverage >= 0.5\n",
    "Coverage >= 0.5 & \nE-value < 1e-3\n",
    "No filter\n",
    "Coverage >= 0.5 & \nE-value < 1e-20\n"),
  Count = c(
    nrow(hits_MinProteins),
    nrow(hits_cov50),
    nrow(hits_cov50_evalue3),
    nrow(hits_any),
    nrow(hits_cov50_evalue20))) %>%
  arrange(desc(Count)) %>%
  mutate(Filter_num = factor(seq_along(`Positive hit filter`)),
         Filter_label = `Positive hit filter`)

# Gradient colors
colors = colorRampPalette(c("lightgreen", "darkgreen"))(nrow(summary_df))
summary_df$color = colors

# Plot
ggplot(summary_df, aes(x = Filter_num, y = Count, fill = Filter_num)) +
  geom_col(show.legend = TRUE, color = "black", size = 0.3, width = 0.7) +
  scale_fill_manual(
    name = "Positive hit filter",
    values = setNames(summary_df$color, summary_df$Filter_num),
    labels = summary_df$Filter_label) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    labels = label_number(scale = 1e-3, suffix = "K")) +
  labs(
    x = "Positive hit filter",
    y = "Number of detected hits in 03/21") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 19),
    legend.text = element_text(size = 16))


################################################################################################
#### Setting the tolerance to gaps cutoffs for all P-P types, 50% coverage, 75% pattern

Table_all_composition_tolerance_to_gaps_cutoff = fread("~/work/PP_profiles_2024/GitHub_structure_0/P-P_detection/Paper_related_data/Supplementary_data/tyPPing_criteria_tables/Table_all_composition_tolerance_to_gaps_cutoff.tsv") 

# Prepare aggregated plot data
plot_data = Table_all_composition_tolerance_to_gaps_cutoff %>%
  mutate(MGE_group = if_else(`MGE type` == `P-P type tested`, "Same type \n(as in subplot title)", "Other type")) %>%
  group_by(`P-P type tested`, `Tolerance to gaps`, MGE_group) %>%
  summarise(`Number of genomes` = sum(`Number of genomes`, na.rm = TRUE), .groups = "drop")

# Horizontal dashed line (max detected P-Ps)
line_data = plot_data %>%
  filter(MGE_group == "Same type \n(as in subplot title)") %>%
  group_by(`P-P type tested`) %>%
  summarise(max_genomes = max(`Number of genomes`, na.rm = TRUE), .groups = "drop")

# Vertical line (tolerance to gaps threshold)
vertical_lines = tibble::tribble(
  ~`P-P type tested`, ~vline_x,
  "AB_1", 7,
  "cp32", 6,
  "N15", 9,
  "P1_1", 11,
  "P1_2", 17,
  "pCAV", 3,
  "pKpn", 2,
  "pMT1", 6,
  "pSLy3", 3,
  "SSU5_pHCM2", 3)

# Plot
ggplot(plot_data, aes(x = `Tolerance to gaps`, y = `Number of genomes`, color = MGE_group)) +
  # Allowed gaps
  geom_rect(
    data = vertical_lines,
    aes(xmin = -Inf, xmax = vline_x, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "green", alpha = 0.12) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  # max detected P-Ps
  geom_hline(
    data = line_data, size = 1,
    mapping = aes(yintercept = max_genomes),
    linetype = "dashed") +
  # tolerance to gaps threshold
  geom_vline(
    data = vertical_lines,
    mapping = aes(xintercept = vline_x),
    color = "red", size = 1) +
  facet_wrap(~ `P-P type tested`, scales = "free_y", nrow = 2) +
  scale_y_continuous(breaks = pretty_breaks(n = 4)) +
  scale_x_continuous(breaks = seq(0, 20, 4), limits = c(-1, 19)) +
  labs(
    x = "Number of allowed missing proteins",
    y = "Number of recovered genomes",
    color = "MGE type") +
  scale_color_manual(values = c("Same type \n(as in subplot title)" = "forestgreen", "Other type" = "tomato")) +
  theme_minimal(base_size = 19) +
  theme(
    strip.text = element_text(face = "bold", size = 19),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 19),
    legend.title = element_text(size = 21),
    legend.text = element_text(size = 19))




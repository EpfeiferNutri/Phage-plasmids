# libraries

library(data.table)
library(scales)
library(readxl)
library(tidyverse)


############### Figure 1A: Profile generation 

# number of genomes used to generate signature profiles
N_genomes =  read_excel( "Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx") %>% 
  filter(`MGE type MM-GRC` %in% c("N15 group", "AB_1 subgroup", "cp32 group", 
                                      "P1_1 subgroup", "P1_2 subgroup", 
                                      "pKpn group", "pSLy3 group", "pCAV group", "pMT1 group", "SSU5_pHCM2 group") & `in 03/21` == "Yes") %>%
  mutate(`P-P type` = str_remove(`MGE type MM-GRC`, " (group|subgroup)$")) %>%
  select(`NCBI accession`, Source, `P-P type`) %>% 
  group_by(Source, `P-P type`) %>% summarise(count = n())

ggplot(N_genomes, aes(fill=Source, y=count, x=`P-P type`)) + 
  geom_bar(position="stack", stat="identity", width = 0.6) +
  coord_flip() + scale_y_reverse() + 
  theme_bw() + 
  ylab("Number of genomes") + xlab("P-P type")+
  theme(strip.background =element_blank(), text = element_text(size = 20)) +
  scale_fill_manual(values = c("#00AAFF", "#FFA64D"))

# number of generated signature profiles
N_signature_genes = read_excel("Publication_related_data/Supplementary_data/T1_Cutoffs_information_table.xlsx")%>% 
  select(`P-P type`, `Total number of signature genes`) 

ggplot(N_signature_genes, aes(y = `Total number of signature genes`, x = `P-P type`)) + 
  geom_bar(position = "stack", stat = "identity", fill = "#38A856", width = 0.6) + 
  scale_y_continuous(breaks = seq(0, 125, 20)) +
  coord_flip() + 
  theme_bw() + 
  xlab("P-P type") + 
  ylab("Number of generated HMMs") + 
  theme(strip.background = element_blank(), text = element_text(size = 20), axis.text.y = element_text(face = "bold"))



############### Figure 1B: Ratio 'Signature genes / All genes'

Ratio_signature_to_all_0321 =  read_excel( "Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx") %>% 
  filter(`MGE type MM-GRC` %in% c("N15 group", "AB_1 subgroup", "cp32 group", 
                                      "P1_1 subgroup", "P1_2 subgroup", 
                                      "pKpn group", "pSLy3 group", "pCAV group", "pMT1 group", "SSU5_pHCM2 group") & `in 03/21` == "Yes") %>%
  mutate(`P-P type` = str_remove(`MGE type MM-GRC`, " (group|subgroup)$")) %>%
  select(`NCBI accession`, Ratio = `Signature to all genes ratio`, `P-P type`)

Ratio_signature_to_all_0321 %>% 
  ggplot(aes(x = `P-P type`, y = Ratio)) +
  geom_boxplot(outlier.shape = NA, fill = "#2ea108", alpha = 0.7) +
  geom_jitter(width = 0.1, color = "darkgreen") +
  theme_bw() + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("P-P type") + 
  ylab("Ratio") +
  theme(strip.background = element_blank(), text = element_text(size = 18)) 


############### Signature_gene_ratio_statistics

Ratio_signature_to_all_0321_stat = Ratio_signature_to_all_0321 %>% 
  group_by(`P-P type`) %>%
  summarise(
  count = n(),
  mean_Ratio = mean(Ratio, na.rm = TRUE),
  median_Ratio = median(Ratio, na.rm = TRUE),
  sd_Ratio = sd(Ratio, na.rm = TRUE),
  min_Ratio = min(Ratio, na.rm = TRUE),
  max_Ratio = max(Ratio, na.rm = TRUE))


############### Figure 1C: P-P score graphs

Profile_information = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx")

ggplot(Profile_information, aes(x=x) ) +
  geom_histogram( aes(x = `P-P score`), bins = 30 , fill="#378B85", col = "#0F595D" ) +
  # geom_label( aes(x=0.3, y=150, label="P-P scores"), color="#378B85", size = 7) +
  geom_histogram( aes(x = `P-P type score`,y=after_stat(-..count..)), bins = 30 , fill= "#88BB00", col = "#4C7500") +
  # geom_label( aes(x=0.3, y=-150, label="P-P type scores"), color="#88BB00", size = 7) +
  theme_bw() + 
  scale_y_continuous(breaks = pretty_breaks(n = 10))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  xlab(" Single profile specificity score values") + ylab("Number of profiles")+
  theme(strip.background =element_blank(), text = element_text(size = 18))


############### Figure 1D: Profile annotation with PHROGs (all types together)

Profile_information = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx") %>% 
  mutate(`PHROG/rep/par annotation` = case_when(!is.na(`PHROGs Annotation`) ~ `PHROGs Annotation`, T ~ `Rep/Par hit`)) %>%
  mutate(Annotation = ifelse(grepl("par", `PHROG/rep/par annotation`, ignore.case = TRUE), "DNA, RNA and nucleotide metabolism: Partition", `Protein category (PHROGs)`)) %>%
  group_by(Annotation) %>% summarise(N = n()) %>% replace_na(list(Annotation = "not annotated")) %>% arrange(N)


ggplot(Profile_information, aes(x = "", y = N, fill = reorder(Annotation, N))) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  labs(
    x = NULL, 
    y = NULL,
    fill = "Protein profile function categories"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    axis.text.y = element_blank()
  ) +
  
  scale_fill_manual(values = c("#6A9EFF", "#66D1C1", "#7FAAFF", "#FFB570", "#D3A8E6", 
                               "#81D6A1", "#7AB8FF", "#FF9A9A", "#70C8FF", "#6DD9B3",  
                               "#D4D4D4"))



############### Figure S1: neff values by P-P type

# NEFF is the diversity of the alignment, calculated as exp of the negative entropy averaged over all columns of the alignment

Profile_information = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx")

Profile_information %>% 
  ggplot(aes(x = factor(`Neff value`))) +  
  geom_bar(col = "darkgreen", fill = "darkgreen", width=0.7) +
  theme_bw() + 
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = c("1.0", "1.1", "1.2", "1.3", "1.4", "1.5")) + 
  facet_wrap(~`P-P type`, scales = "free_y", nrow = 2) +
  xlab("Neff value by P-P type") + 
  ylab("Number of profiles") +
  theme(strip.background = element_blank(), text = element_text(size = 17))



############### Figure S2: Profile scores by P-P type 

Profile_information = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx")

ggplot(Profile_information, aes(x=x) ) +
  geom_histogram( aes(x = `P-P score`), bins = 30 , fill="#378B85", col = "#0F595D" ) +
  # geom_label( aes(x=0.3, y=150, label="P-P scores"), color="#378B85", size = 7) +
  geom_histogram( aes(x = `P-P type score`,y=after_stat(-..count..)), bins = 30 , fill= "#88BB00", col = "#4C7500") +
  # geom_label( aes(x=0.3, y=-150, label="P-P type scores"), color="#88BB00", size = 7) +
  theme_bw() + 
  facet_wrap(~`P-P type`, scales="free_y", nrow = 2) +
  scale_y_continuous(breaks = pretty_breaks(n = 6))+
  scale_x_continuous(breaks = seq(0,1,0.3))+
  xlab(" Single profile specificity score values") + ylab("Number of profiles")+
  theme(strip.background =element_blank(), text = element_text(size = 18))



############### Figure S3: Functional annotation of signature profiles by P-P type (with PHROGs)

Profile_information = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx")

Profile_information_stat_by_type = Profile_information %>%
  mutate(`PHROG/rep/par annotation` = case_when(!is.na(`PHROGs Annotation`) ~ `PHROGs Annotation`, T ~ `Rep/Par hit`)) %>%
  mutate(Annotation = ifelse(grepl("par", `PHROG/rep/par annotation`, ignore.case = TRUE), "DNA, RNA and nucleotide metabolism: Partition", `Protein category (PHROGs)`)) %>%
  group_by(Annotation, `P-P type`) %>% summarise(N = n()) %>% replace_na(list(Annotation = "not annotated")) %>% arrange(N)

ggplot(Profile_information_stat_by_type, aes(x = "", y = N, fill = reorder(Annotation, N))) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  labs(
    x = NULL, 
    y = NULL,
    fill = "Protein profile function categories"
  ) +
  facet_wrap(~`P-P type`, nrow = 2) +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    axis.text.y = element_blank()
  ) +
  
  scale_fill_manual(values = c("#6A9EFF", "#66D1C1", "#7FAAFF", "#FFB570", "#D3A8E6", 
                               "#81D6A1", "#7AB8FF", "#FF9A9A", "#70C8FF", "#6DD9B3",  
                               "#D4D4D4"))


########## Tables


All_MGE_information_table = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

Profile_information_table = read_excel("Publication_related_data/Supplementary_data/S2_Profile_information_table.xlsx")


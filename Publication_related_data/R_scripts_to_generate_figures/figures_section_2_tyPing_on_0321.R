# libraries

library(data.table)
library(scales)
library(readxl)
library(tidyverse)


############### 03/21 validation data

All_MGE_information_table = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

Cutoffs_information = read_excel("Publication_related_data/Supplementary_data/T1_Cutoffs_information_table.xlsx")

FINAL_PREDICTION_ = fread("Publication_related_data/Supplementary_data/Final_prediction_table_0321_by_tyPPing.tsv") 

# not used
# FINAL_SUMMARY = fread("Publication_related_data/Supplementary_data/All_hmm_hits_table_0321_by_tyPPing.tsv")

FINAL_PREDICTION = FINAL_PREDICTION_ %>% 
  full_join(All_MGE_information_table %>% filter(`in 03/21` == "Yes"), by=c("Genome ID"="NCBI accession")) %>%
  left_join(Cutoffs_information %>% select(`P-P type`, `Lower size limit (bp)`, `Upper size limit (bp)`)) %>% group_by(`Genome ID`)


################# Figure 2B: Venn diagram numbers

# MinProteins U Composition
aUb = FINAL_PREDICTION %>% filter(!is.na(`Predicted by`)) %>% summarise(n=n())

# MinProteins
a = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition" | `Predicted by` == "MinProteins") %>% summarise(n=n())

# Composition
b = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition" | `Predicted by` == "Composition") %>% summarise(n=n())

# MinProteins & Composition 
aQb = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition") %>% summarise(n=n())

# MinProteins & Size
aQsize = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition" | `Predicted by` == "MinProteins") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())

# Composition & Size
bQsize = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition" | `Predicted by` == "Composition") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())

# MinProteins & Composition & Size
aQbQsize = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())

# MinProteins & Composition & !Size
aQb_size = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition") %>%
  filter (`Genome size (bp)` < `Lower size limit (bp)` | `Genome size (bp)` > `Upper size limit (bp)`) %>% summarise(n=n())

# MinProteins only & !Size
a_ = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins") %>% 
  anti_join(FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition"), by ="Genome ID") %>%
  filter (`Genome size (bp)` < `Lower size limit (bp)` | `Genome size (bp)` > `Upper size limit (bp)`) %>% summarise(n=n())

# MinProteins only & Size
a_size = FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins") %>% 
  anti_join(FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition"), by ="Genome ID") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())

# Composition only & !Size
b_ = FINAL_PREDICTION %>% filter(`Predicted by` == "Composition") %>% 
  anti_join(FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition"), by ="Genome ID") %>%
  filter (`Genome size (bp)` < `Lower size limit (bp)` | `Genome size (bp)` > `Upper size limit (bp)`) %>% summarise(n=n())

# Composition only & Size
b_size = FINAL_PREDICTION %>% filter(`Predicted by` == "Composition") %>% 
  anti_join(FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition"), by ="Genome ID") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())


################# Figure 2C: Confidence levels

# Number of confident predictions 
FINAL_PREDICTION_conf = FINAL_PREDICTION %>% filter(`Confidence level` != "Low")

FINAL_PREDICTION_conf_high = FINAL_PREDICTION_conf %>% filter(`Confidence level` == "High")
FINAL_PREDICTION_conf_medium = FINAL_PREDICTION_conf %>% filter(`Confidence level` == "Medium")

# Number of low-confident predictions
FINAL_PREDICTION_conf_low = FINAL_PREDICTION %>% filter(`Confidence level` == "Low")
FINAL_PREDICTION_conf_low_short  = FINAL_PREDICTION_conf_low %>% filter(`Genome size (bp)` < `Lower size limit (bp)`)

Low_conf_summary = FINAL_PREDICTION_conf_low %>% 
  group_by(`P-P type`) %>% summarise(n=n())

######################## Figure 2C: Number of P-P with confidence levels 

confidence_levels_count = data.frame(
  Confidence_level = c("High", "Medium","Low"), 
  Count = c(nrow(FINAL_PREDICTION_conf_high), nrow(FINAL_PREDICTION_conf_medium), nrow(FINAL_PREDICTION_conf_low)))

# Set factor levels to maintain the desired order
confidence_levels_count$Confidence_level <- factor(confidence_levels_count$Confidence_level, levels = c("High", "Medium","Low"))
confidence_colors = c("High" = "#51ad13", "Medium" = "#257ec2", "Low" = "#d4522f")


ggplot(confidence_levels_count, aes(x = Confidence_level, y = Count, fill = Confidence_level)) +
  geom_bar(stat = "identity", width = 0.7) +
  # coord_flip() +  # Horizontal bars
  scale_fill_manual(values = confidence_colors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),  
    axis.ticks = element_blank(), 
    panel.grid = element_blank(),  
    text = element_text(size = 20)
  )


################# Figure 2C: Comparison to MM-GRC

overlap = FINAL_PREDICTION_conf %>% filter(!`MGE type MM-GRC` %in% c("Plasmid", "Phage"))

tyPPing_only = FINAL_PREDICTION_conf %>% filter(`MGE type MM-GRC` %in% c("Plasmid", "Phage"))

MM_GRC_only = anti_join(FINAL_PREDICTION %>% filter(`MGE type MM-GRC` %in% c("N15 group", "AB_1 subgroup", "cp32 group", 
                                                                               "P1_1 subgroup", "P1_2 subgroup", "P11 group","P12 group", 
                                                                               "pKpn group", "pSLy3 group", "pCAV group", "pMT1 group", "SSU5_pHCM2 group") ), FINAL_PREDICTION_conf, by="Genome ID")



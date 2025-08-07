# libraries

library(data.table)
library(scales)
library(readxl)
library(tidyverse)


############### 0523 test data

All_MGE_information_table = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

Cutoffs_information = read_excel("Publication_related_data/Supplementary_data/T1_Cutoffs_information_table.xlsx")

FINAL_PREDICTION_ = fread("Publication_related_data/Supplementary_data/Final_prediction_table_0523_by_tyPPing.tsv")
# not used
# FINAL_SUMMARY = fread("Publication_related_data/Supplementary_data/All_hmm_hits_table_0523_by_tyPPing.tsv")

FINAL_PREDICTION = FINAL_PREDICTION_ %>%
  full_join(All_MGE_information_table, by=c("Genome ID"="NCBI accession"))  %>% filter(`in 03/21` == "No") %>%
  left_join(Cutoffs_information %>% select(`P-P type`, `Lower size limit (bp)`, `Upper size limit (bp)`)) %>% group_by(`Genome ID`)


################# Venn diagram numbers

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

# Composition only & !Size
b_size = FINAL_PREDICTION %>% filter(`Predicted by` == "Composition") %>% 
  anti_join(FINAL_PREDICTION %>% filter(`Predicted by` == "MinProteins + Composition"), by ="Genome ID") %>%
  filter (`Genome size (bp)` >= `Lower size limit (bp)` & `Genome size (bp)` <= `Upper size limit (bp)`) %>% summarise(n=n())


################# Figure 2E: Confidence levels

# Number of confident predictions 
FINAL_PREDICTION_conf = FINAL_PREDICTION %>% filter(`Confidence level` != "Low")

FINAL_PREDICTION_conf_high = FINAL_PREDICTION_conf %>% filter(`Confidence level` == "High")
FINAL_PREDICTION_conf_medium = FINAL_PREDICTION_conf %>% filter(`Confidence level` == "Medium")

# Number of low-confident predictions
FINAL_PREDICTION_conf_low = FINAL_PREDICTION %>% filter(`Confidence level` == "Low")
FINAL_PREDICTION_conf_low_short  = FINAL_PREDICTION_conf_low %>% filter(`Genome size (bp)` < `Lower size limit (bp)`)

Low_conf_summary = FINAL_PREDICTION_conf_low %>% 
  group_by(`P-P type`) %>% summarise(n=n())


######################## Figure 2E: Number of P-Ps with confidence levels 

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


################# Figure 3E: comparison to MM-GRC

overlap = FINAL_PREDICTION_conf %>% filter(!`MGE type MM-GRC` %in% c("Plasmid", "Phage"))

tyPPing_only = FINAL_PREDICTION_conf %>% filter(`MGE type MM-GRC` %in% c("Plasmid", "Phage"))

MM_GRC_only = anti_join(FINAL_PREDICTION %>% filter(`MGE type MM-GRC` %in% c("N15 group", "AB_1 subgroup", "cp32 group", 
                                                                                     "P1_1 subgroup", "P1_2 subgroup", "P11 group","P12 group", 
                                                                                     "pKpn group", "pSLy3 group", "pCAV group", "pMT1 group", "SSU5_pHCM2 group") ), FINAL_PREDICTION_conf, by="Genome ID")


################# Figure 2D: Number of typed P-Ps in 05/23 in comparison to 03/21

# known P-Ps from 0321 detected by tyPPing
FINAL_PREDICTION_0321 = fread("Publication_related_data/Supplementary_data/Final_prediction_table_0321_by_tyPPing.tsv") %>%
  filter(`Confidence level` != "Low") %>% select(`P-P type`) %>% 
  mutate(`P-P type` = case_when(`P-P type` == "SSU5_pHCM2" ~ "SSU5", T ~ `P-P type`)) %>% 
  group_by(`P-P type`) %>% summarise("03/21"=n())


detected_PPs = FINAL_PREDICTION_conf %>% select(`P-P type`) %>% 
  mutate(`P-P type` = case_when(`P-P type` == "SSU5_pHCM2" ~ "SSU5", T ~ `P-P type`)) %>% 
  group_by(`P-P type`) %>% summarise("05/23"=n()) %>%
  left_join(FINAL_PREDICTION_0321) %>% mutate("0321" = ifelse(`P-P type` == "P1_1", `03/21`-2, `03/21`))

# Convert to long format for stacking
detected_PPs_ = detected_PPs %>% pivot_longer(cols = c("03/21", "05/23"), names_to = "Database", values_to = "Number of detected P-Ps")

detected_PPs_$Database = factor(detected_PPs_$Database, levels = c("05/23", "03/21"))

custom_colors = c("03/21" = "grey50", "05/23" = "#81b55c")  

ggplot(detected_PPs_, aes(x = `P-P type`, y = `Number of detected P-Ps`, fill = Database)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_bw() +
  coord_flip() +  # Horizontal bars
  labs(x = "P-P Type", y = "Number of detected P-Ps") +
  theme(text = element_text(size = 24), axis.text.y = element_text(face = "bold"))  






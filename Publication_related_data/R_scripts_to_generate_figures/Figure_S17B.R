library(data.table)
library(readxl)
library(UpSetR)
library(tidyverse)


######################################### test different vConTACT2 modes

genome_by_genome_overview_0523 = fread('Publication_related_data/R_scripts_to_generate_figures/vConTACT2_genome_by_genome_overview_for_0523_plasmids.csv')

# vConTACT2 default clustering

reference_0321_PP_list_for_vConTACT2 = read_excel("geNomad_vConTACTv2/vConTACTv2/reference_0321_PP_list_for_vConTACT2.xlsx") 

PPs_only_Summary_0523_vcontact = genome_by_genome_overview_0523 %>% 
  left_join(reference_0321_PP_list_for_vConTACT2, by=c("Genome"="NCBI accession")) %>%
  select(Genome, `VC Status`, VC, `VC Subcluster`,`P-P type`)

PPs_only_Summary_0523_vcontact_1 = PPs_only_Summary_0523_vcontact %>% 
  # `Virus cluster` is the most general cluster 
  mutate(`Virus cluster` = case_when(`VC Status` == "Singleton" ~ str_c(row_number()," Singleton"), 
                                     `VC Status` == "Outlier" ~ str_c(row_number()," Outlier"),
                                     `VC Status` == "Clustered" ~ str_remove(`VC`, pattern="_[0-9]+$"), 
                                     str_detect(`VC Status`, "Overlap") ~ str_remove_all(`VC Status`, pattern = "^Overlap \\(|\\)$"), 
                                     `VC Status` == "Clustered/Singleton" ~ str_remove(`VC`, pattern="_[0-9]+$"),
                                     T ~ NA),
         `P-P type vConTACT2` = case_when(!is.na(`P-P type`) ~ `P-P type`, 
                                          T ~ NA)) %>%
  # prediction by vcontact;  if conflicted groups in one cluster --> all are kept
  group_by(`Virus cluster`) %>%
  mutate(`P-P type vConTACT2` = if (all(is.na(`P-P type`))) NA else paste(unique(na.omit(`P-P type`)), collapse = "; ")) %>%
  ungroup() %>%
  filter(!is.na(`P-P type vConTACT2`)) 


# vConTACT2 clustered only

reference_0321_PP_list_for_vConTACT2 = read_excel("geNomad_vConTACTv2/vConTACTv2/reference_0321_PP_list_for_vConTACT2.xlsx")

PPs_only_clustered_Summary_0523_vcontact = genome_by_genome_overview_0523 %>% 
  filter(`VC Status` == "Clustered") %>%
  left_join(reference_0321_PP_list_for_vConTACT2, by=c("Genome"="NCBI accession")) %>%
  select(Genome, `VC Status`, VC, `VC Subcluster`,`P-P type`)

PPs_only_clustered_Summary_0523_vcontact_3 = PPs_only_clustered_Summary_0523_vcontact %>% 
  # `Virus cluster` is the most general cluster 
  mutate(`Virus cluster` =  str_remove(`VC`, pattern="_[0-9]+$"),
         `P-P type vConTACT2` = case_when(!is.na(`P-P type`) ~ `P-P type`, 
                                          T ~ NA)) %>%
  # prediction by vcontact;  if conflicted groups in one cluster --> all are kept
  group_by(`Virus cluster`) %>%
  mutate(`P-P type vConTACT2` = if (all(is.na(`P-P type`))) NA else paste(unique(na.omit(`P-P type`)), collapse = "; ")) %>%
  ungroup() %>%
  filter(!is.na(`P-P type vConTACT2`)) 


####################################################################################################


all_mge = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx") %>% 
  filter(`in 05/23` == "Yes", `in 03/21` == "No") %>%
  select(`NCBI accession`, `Replicon size (bp)`,  `MGE type MM-GRC`)

ALL_mge = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx")

all_test_vcontact = all_mge %>% 
  left_join(PPs_only_Summary_0523_vcontact_1 %>% select(`NCBI accession` = Genome, `vConTACT2 PPs_only` = `P-P type vConTACT2`)) %>%
  left_join(PPs_only_clustered_Summary_0523_vcontact_3 %>% select(`NCBI accession` = Genome, `vConTACT2 PPs_only_clustered` = `P-P type vConTACT2`)) 

all_test_vcontact_plot_ = all_test_vcontact %>%
  mutate(`MM-GRC` = case_when(`MGE type MM-GRC` != "Plasmid" & `MGE type MM-GRC` != "Phage" ~ 1, T ~ 0), 
         `vConTACT2 default` = case_when(!is.na(`vConTACT2 PPs_only`) ~ 1, T ~ 0), 
         `vConTACT2 clustered` = case_when(!is.na(`vConTACT2 PPs_only_clustered`) ~ 1, T ~ 0))

all_test_vcontact_plot = all_test_vcontact_plot_ %>%
  select(`NCBI accession`, `MM-GRC`, `vConTACT2 default`, `vConTACT2 clustered`) %>% 
  column_to_rownames("NCBI accession")

upset(data = all_test_vcontact_plot,
      nsets = 6,
      mainbar.y.label = "Putative P-Ps",
      order.by = "freq",
      sets.x.label = "", 
      decreasing = TRUE, 
      main.bar.color = "skyblue",
      sets.bar.color = "orange",
      text.scale = 1.6,
      point.size = 4)

colSums(all_test_vcontact_plot)

binary_cols = c("MM-GRC", "vConTACT2 default", "vConTACT2 clustered")

all_test_vcontact_plot_$combo_id = apply(all_test_vcontact_plot_[, binary_cols], 1, paste0, collapse = "_")

combination_dfs = split(all_test_vcontact_plot_, all_test_vcontact_plot_$combo_id)

subset_1 = combination_dfs[["1_1_1"]] %>% select(`NCBI accession`) %>% left_join(ALL_mge)
subset_2 = combination_dfs[["0_1_1"]] %>% select(`NCBI accession`) %>% left_join(ALL_mge)
subset_3 = combination_dfs[["0_1_0"]] %>% select(`NCBI accession`) %>% left_join(ALL_mge)
subset_4 = combination_dfs[["1_0_0"]] %>% select(`NCBI accession`) %>% left_join(ALL_mge)
subset_5 = combination_dfs[["1_1_0"]] %>% select(`NCBI accession`) %>% left_join(ALL_mge)

######################################

# libraries
library(tidyverse)
library(readxl)

###
# Path to table S1, with all MGEs and P-P predictions
s1_all_mge = read_excel("Table_S1_All_MGEs.xlsx")

plasmids_0523_genomad = s1_all_mge %>% filter(`MGE type geNomad 0523` == "Plasmid")
phage_0523_genomad = s1_all_mge %>% filter(`MGE type geNomad 0523` == "Phage")

s1_all_mge_signature = s1_all_mge %>% filter(!is.na(`Signature to all genes ratio`)) %>% 
  dplyr::count(`in 0321`, `in 0523`, `MGE type (by Model)`)

s1_all_mge_low_0321 = s1_all_mge %>% filter(`Confidence level (tyPPing)` == "Low", `in 0321` == "Yes")

total_pps = s1_all_mge %>% 
  # take only the detectd ones
  filter((Source =="Plasmid_DB"&`MGE type geNomad 0523` == "Phage")  | !(`MGE type MM-GRC` %in% c("Phage","Plasmid")) | `Confidence level (tyPPing)` %in% c("High","Medium") ) %>%
  mutate(genomad_plasmid = ifelse(`NCBI accession` %in% plasmids_0523_genomad$`NCBI accession`,1,0 ))
  # only on the group level

total_pps_grouped= total_pps %>%
  filter(grepl(`MGE type MM-GRC`, pattern = "group") | grepl(`P-P type tyPPing`, pattern = "group")|grepl(`P-P type 0523 geNomad+vConTACT2`, pattern = "group"))

##singletons by MM-GRC
singletons_MMGRC = total_pps %>% group_by(`MGE type MM-GRC`) %>% mutate(n_members_MMGRC = n()) %>%
  filter(`in 05/23` == "Yes", `in 03/21` == "No", n_members_MMGRC == 1) %>% ungroup() 

# singletons by geNomad
genomad = total_pps %>% select(Accession = `NCBI accession`, genomad_MGE = `MGE type geNomad 0523`, Source, `in 03/21`, `in 05/23`) %>%
  mutate(genomad_phage = ifelse( Source == "Plasmid_DB"& genomad_MGE=="Phage",1,0), genomad_prophage = ifelse(grepl(genomad_MGE,pattern="Prophage"),1,0) 
         ,genomad_plasmid = ifelse(genomad_MGE == "Plasmid",1,0))

#
vcontacts = total_pps %>% 
  filter( Source == "Plasmid_DB"& `MGE type geNomad 0523`=="Phage") %>%
  select(Accession = `NCBI accession`, vcontact_grouping = `P-P type 0523 geNomad+vConTACT2`) %>%
  mutate(vcontact = ifelse(grepl(vcontact_grouping, pattern = "group|community|Singleton"), 1, 0)) %>%
  filter(vcontact_grouping != "NA")

#
singletons_genomad = genomad %>% left_join(vcontacts) %>%
  filter(`in 05/23` == "Yes", `in 03/21` == "No" ) %>%
  replace_na( replace = list("vcontact"=0)) %>%
  filter(genomad_phage ==1, vcontact == 0)

### put singletons together
singletons_upset_df = 
  full_join(singletons_genomad, singletons_MMGRC %>% select("Accession"=`NCBI accession`,mmgrc=n_members_MMGRC), by = "Accession") %>%
  select(Accession, genomad_phage, mmgrc) %>% replace(is.na(.),0) %>% column_to_rownames("Accession")

# plot (Figure 3C)
  UpSetR::upset(#data=df_to_plot
      data=singletons_upset_df
      #,nsets = 10
      ,mainbar.y.label = "putative P-P singletons"
      #,show.numbers = T
      ,order.by = "freq"
      ,sets.x.label = "
      " 
      #,decreasing = T 
      ,main.bar.color = alpha("tomato", alpha = 0.8),
      sets.bar.color = "orange",
      text.scale = 1.8,
      point.size = 6)


## PPs by MMGRC  
MMGRC = total_pps %>% select(Accession = `NCBI accession`, mge_type_model = `MGE type MM-GRC`) %>% 
  mutate(mmgrc = ifelse(mge_type_model != "Plasmid", 1, 0) ) 

# PPs by tyPPing 
tyPPings = total_pps %>% select(Accession = `NCBI accession`, tyPPing_confidence=`Confidence level (tyPPing)`) %>%
  mutate(tyPPing = ifelse(tyPPing_confidence %in% c("High","Medium"), 1, 0))

### Merge table and make upset plot

## put together and make upset plot
pre_upset_plot = tyPPings %>% left_join(vcontacts) %>% 
    left_join(MMGRC) %>% left_join(genomad) %>% replace_na( replace = list("vcontact"=0)) %>%
    group_by(mge_type_model) %>% mutate(n_members_MMGRC = n()) %>% ungroup()
  
#  
upset_plot_df = pre_upset_plot %>%
  filter(`in 05/23` == "Yes", `in 03/21` == "No" ) %>%
    select(Accession, mmgrc,geno_vcontact=vcontact,tyPPing) %>% 
    anti_join(singletons_upset_df %>% rownames_to_column("Accession"), by ="Accession") %>%
    column_to_rownames("Accession") 

# plot (Figure 3B)
  UpSetR::upset(#data=df_to_plot
      data=upset_plot_df
      ,nsets = 10
      ,mainbar.y.label = "putative P-Ps"
      #,show.numbers = T
      ,order.by = "freq"
      ,sets.x.label = "
      " 
      #,decreasing = T 
      ,main.bar.color = alpha("tomato", alpha = 0.8),
      sets.bar.color = "orange",
      text.scale = 1.8,
      point.size = 6)
  
### 
  
  ### Number of P-Ps classed by geNomad as plasmids and intrgrated prophages
  
  # filter for tbl
  tyPPing_geno = pre_upset_plot %>% filter(`in 05/23` == "Yes", `in 03/21` == "No", tyPPing_confidence %in% c("High","Medium") ) %>% dplyr::count(genomad_MGE,Source) %>% mutate(Method = "tyPPing")
  mmgrc_geno = pre_upset_plot %>% filter(`in 05/23` == "Yes", `in 03/21` == "No", mmgrc == 1 ) %>% dplyr::count(genomad_MGE,Source) %>% mutate(Method = "MM-GRC")
  
  # make the count tbl
  bar_plot_df = bind_rows(tyPPing_geno,mmgrc_geno) %>% 
    mutate(geNomad_predict = ifelse(grepl(genomad_MGE, pattern="Plasmid|Prophage"), "Plasmid or Prophage",genomad_MGE)) %>%
    group_by(Method,geNomad_predict) %>% summarise(n=sum(n)) %>%
    group_by(Method) %>% mutate(Fraction = round(100*n/sum(n),1)) %>% ungroup()
  
  # barplot, # Fig S17 A
  ggplot(bar_plot_df, aes( x = n, y = Method, fill = geNomad_predict)) +
  geom_col() +
  scale_fill_manual(values = c("tomato","lightblue","grey")) +
  theme_bw()+theme(text = element_text(size = 14), legend.position = "top", legend.title = element_blank())




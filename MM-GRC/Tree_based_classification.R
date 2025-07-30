### required packages ###

### Uncomment and install if not installed:
# library(devtools)
# install_github('arendsee/rhmmer')

library("ranger")
library("readxl")
library("rhmmer")
library("data.table")
library("pROC")
library(dtplyr)
library("foreach")
library("tidyverse")
library("seqinr")


# CHANGE ALL THE PATHS below to match your local file system

### ------------------
### DATA PREPROCESSING 
### ------------------

### 1. Load functional annotation tables ####

# annotation category 1: pVOGs with TIGRFAM/PFAM
func_cat1 = read_excel(path = "path_to/P-P_detection/MM-GRC/200122_func_cat_pvogs_pfams_tigfams_v4.xlsx", sheet = 1) %>% 
  select(HMM, Functional_categorization)
# annotation category 2: abundant pVOGs that cluster and do not cluster with each other (VQ: 0.75, proteins: 15)
func_cat2 = read_excel(path = "path_to/P-P_detection/MM-GRC/200122_func_cat_pvogs_pfams_tigfams_v4.xlsx", sheet = 2) %>% 
  select(HMM, Functional_categorization_2)

# Merge both into a unified table
func_cat = full_join(func_cat1,func_cat2)


### 2. Process HMMER --domtblout output

# Read HMMER search output (domtblout format)
HMM_output_df = data.frame(NULL)
tmp_df = rhmmer::read_domtblout("path_to_hmm_search_output_for_plasmid_proteins_against_phage_profiles/hmm_output.tbl.out") 

# Filter hits: domain e-value < 1e-3 and domain alignment coverage > 50%
tmp_df = tmp_df %>% filter(domain_ievalue < 1e-3) %>%
  mutate( ali_diff = ali_to+1 - ali_from, cov_ali_prot = ali_diff/domain_len, cov_ali_hmm = ali_diff/qlen) %>%
  filter(cov_ali_hmm > 0.5)

HMM_output_df = HMM_output_df %>% bind_rows(tmp_df) 

# Extract HMM ID and DB source (PFAM, TIGRFAM, or pVOG)
HMM_output_df = HMM_output_df %>% mutate(
  HMM = case_when(!is.na(query_accession) ~ query_accession, TRUE ~ query_name),
  HMM = str_split(HMM,pattern="\\.", simplify = T)[,1],
  HMM_DB = case_when(grepl(HMM, pattern = "VOG")~"pVOG",grepl(HMM, pattern = "TIGR")~"TIGFRFAM",grepl(HMM, pattern = "PF")~"PFAM")) 

HMM_output_df_ = HMM_output_df %>% select(protein_id=domain_name, HMM,HMM_DB, cov_ali_hmm, domain_ievalue)


### OPTIONAL: if the dataset is too big to process ### 
# 1. save processed HMM output
write_tsv(HMM_output_df_, "path_to_HMM_output_df_/filtered_pfam_tigr_pvog_plasmids_0523.txt")
# 2. read processed HMM output (new session)
HMM_output_df_ = fread("path_to_HMM_output_df_/filtered_pfam_tigr_pvog_plasmids_0523.txt")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Keep only HMMs present in functional annotation table
HMM_output_df = HMM_output_df_  %>%  semi_join(func_cat, by = "HMM") 


### 3. Load protein sequences (FASTA)
# protein sequences should be named as given by default from prodigal: 
# genome/contig-name_protein-number
# e.g. contig1_1, genome/contig1_2, genome/contig1_3, etc.

plasmids_prots_lst = read.fasta("path_to_plasmid_protein_sequences_fasta_file/plasmids.prt", seqtype = "AA", as.string = T)

# Extract protein metadata
plasmids_prots_df = data.frame( protein_id = unlist(lapply(plasmids_prots_lst,getName)),
                                       protein_size = unlist(lapply(plasmids_prots_lst,getLength))) 

plasmids_prots_df = plasmids_prots_df %>%
  mutate(protein_id = str_split(protein_id, pattern = "\t", simplify = T)[,1],
         replicon = str_remove(protein_id, pattern = "_[0-9]{1,5}$"))

rm(plasmids_prots_lst)
gc()

# Count proteins per replicon
plasmid_prt_nr = plasmids_prots_df %>% group_by(replicon) %>% summarise(pr_nr = n())

### 4. Merge HMM hits with protein info

HMM_prots_hits_plasmids = plasmids_prots_df %>% left_join(HMM_output_df) %>% 
  replace(is.na(.), "not_hit") 

# Summary table
table(HMM_prots_hits_plasmids$HMM_DB)


### -----------------------------
###  FEATURE EXTRACTION & 
###  RF MODEL PREDICTION PIPELINE
### -----------------------------

### 1. Generate feature table

# Join HMM hit data with functional categories and summarize by replicon
Model_features_df = HMM_prots_hits_plasmids %>% left_join(func_cat) %>% 
  replace(., is.na(.), "N_A") %>% group_by(replicon) %>% 
  summarize(  
    DNA_met_Reg_Rec1 = sum(grepl("DNA",Functional_categorization  ))
    ,Lysis1=sum(grepl("Lysis",Functional_categorization))                   
    ,Structure1 = sum(grepl("Structure", Functional_categorization  ))
    ,Pack_Ass_Inj1=sum(grepl("Pack",Functional_categorization ))                                                       
    ,Others1 = sum(grepl("Others",    Functional_categorization  ))
    ,Unknown1 = sum(grepl("Unknown", Functional_categorization ))
    ,DNA_met_Reg_Rec2 = sum(grepl("DNA", Functional_categorization_2  ))
    ,Lysis2 =sum(grepl("Lysis", Functional_categorization_2))                   
    ,Structure2 =sum(grepl("Structure", Functional_categorization_2  ))
    ,Pack_Ass_Inj2 =sum(grepl("Pack", Functional_categorization_2))                                                       
    ,Others2 =sum(grepl("Others",    Functional_categorization_2  ))
    ,Unknown2 = sum(grepl("Unknown", Functional_categorization_2)) 
    ,DNA_met_Reg_Rec3 = length(grep("DNA",       c(Functional_categorization,Functional_categorization_2)))
    ,Lysis3           = length(grep("Lysis",     c(Functional_categorization,Functional_categorization_2)))                   
    ,Structure3       = length(grep("Structure", c(Functional_categorization,Functional_categorization_2)))
    ,Pack_Ass_Inj3    = length(grep("Pack",      c(Functional_categorization,Functional_categorization_2)))                                                       
    ,Others3          = length(grep("Others",    c(Functional_categorization,Functional_categorization_2)))
    ,Unknown3         = length(grep("Unknown",   c(Functional_categorization,Functional_categorization_2))) 
    ,prots = length(unique(protein_id))
    ,HMM_hits = sum(!grepl(HMM,pattern="not_hit"))
    ,VOG_hits = sum(grepl(HMM,pattern="VOG"))
    ,PFAM_hits = sum(grepl(HMM,pattern="PF"))
    ,TIGR_hits = sum(grepl(HMM,pattern="TIGR"))
    ,VOGcat1_hits = sum(!grepl(Functional_categorization,pattern="N_A"))
    ,VOGcat2_hits = sum(!grepl(Functional_categorization_2,pattern="N_A"))
    ,unhit = sum(grepl(HMM,pattern="not_hit"))
    ,phage_hmm = VOGcat1_hits+ VOGcat2_hits+ TIGR_hits + PFAM_hits
    #,Hit_ratio = HMM_hits/prots
    ,unhit_ratio = unhit/prots
    ,VOG_ratio = VOG_hits/prots
    ,PFAM_ratio = PFAM_hits/prots
    ,TIGR_ratio = TIGR_hits/prots
    #Cat1
    ,DNA_met_ratio1 = DNA_met_Reg_Rec1/prots
    ,Lysis_ratio1 = Lysis1/prots
    ,Structure_ratio1 = Structure1/prots
    ,Pack_Ass_ratio1 = Pack_Ass_Inj1/prots
    ,Others_ratio1 = Others1/prots
    ,Unknown_ratio1 = Unknown1/prots
    #Cat2
    ,DNA_met_ratio2 = DNA_met_Reg_Rec2/prots
    ,Lysis_ratio2 = Lysis2/prots
    ,Structure_ratio2 = Structure2/prots
    ,Pack_Ass_ratio2 = Pack_Ass_Inj2/prots
    ,Others_ratio2 = Others2/prots
    ,Unknown_ratio2 = Unknown2/prots
    ,Cat1_ratio = VOGcat1_hits/prots
    ,Cat2_ratio = VOGcat2_hits/prots
  ) %>% ungroup()


### OPTIONAL: if the dataset is too big to process ### 
# 1. save feature table
write_tsv(Model_features_df, "path_to_Model_features_df/model_features_unpredicted.txt")
# 2. read feature table (new session)
Model_features_df = fread("path_to_Model_features_df/model_features_unpredicted.txt")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### 2. Apply the pretrained RF models ###

# Read in the models
result_df = Model_features_df %>% select(replicon) 
models = dir("path_to/P-P_detection/MM-GRC/models/",pattern="pp_model")

# Predict phage/plasmid probability with each model
foreach(i=seq_along(models)) %do% {
  phage_plasmid_model = readRDS(file = str_c("path_to_Model_features_df/MM-GRC/models/",models[i]))
  # make the prediction
  tmp_model_features = Model_features_df 
  tmp_model_features$pred_stat = predict(object = phage_plasmid_model, data = Model_features_df)$predictions 
  tmp_model_features = tmp_model_features %>% dplyr::select(replicon, pred_stat)
  result_df = result_df %>% left_join(tmp_model_features, by= "replicon")
}

colnames(result_df) = c("replicon", str_c(c("psc"), sort(rep(1:10,1))))

### -----------------------------
###  SUMMARIZE PREDICIONS &
###  SAVE RESULTS
### -----------------------------

### 1. Compute final RF model scores

stats_df = result_df %>% mutate(mean_prob=apply(select(.,starts_with("psc")),1,mean), sd_prob = apply(select(., starts_with("psc")),1,sd)) %>% 
  select(replicon,mean_prob,sd_prob)

### Replicons with >0.5 probability to be P-Ps
### Note: only replicons with size <= 300kb & size >= 10kb should be considered 
predicted_PP_list = stats_df %>% filter(mean_prob > 0.5)

# save the prediction results
write_tsv(stats_df, "path_to_model_stats/model_features_predicted_stats.txt")
write_tsv(predicted_PP_list, "path_to_predicted_PP_list/PP_list_stats_by_model.tsv")


### 2. Extract and save protein sequences (for P-Ps)

# Read and process protein FASTA file
plasmids_prots_lst = read.fasta("path_to_plasmid_protein_sequences_fasta_file/plasmids.prt", seqtype = "AA", as.string = T)

protein_seqs = data.frame(  protein_id = unlist(lapply(plasmids_prots_lst,getName))
                            ,protein_size = unlist(lapply(plasmids_prots_lst,getLength))
                            ,header = unlist(lapply(plasmids_prots_lst,getAnnot))
                            ,protein_seq = unlist(lapply(plasmids_prots_lst,getSequence, as.string= T))) %>%         
  mutate(protein_id = str_split(protein_id, pattern = "\t", simplify = T)[,1],
                     replicon = str_remove(protein_id, pattern = "_[0-9]{1,5}$"))


# Filter for predicted P-P replicons
pp_protein_seqs = protein_seqs %>% semi_join(predicted_PP_list, by = "replicon")

# Save protein sequences as FASTA
write.fasta(sequences = as.list(pp_protein_seqs$protein_seq), names =  as.list(str_remove(pp_protein_seqs$header, pattern = "^>"))
            ,file.out = "path_to_predicted_PP_protein_sequences_fasta_file/PP_protein_seqs.fasta")




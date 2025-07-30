### required packages ###

library(data.table)
library(seqinr)
library(tidyverse)

###   1. Load protein sequences (FASTA)
# protein sequences should be named as given by default from prodigal: 
# genome/contig-name_protein-number
# e.g. contig1_1, genome/contig1_2, genome/contig1_3, etc.

prot_info_lst = read.fasta(seqtype = "AA", as.string = T, file="pp_prots.fasta")

# Extract protein IDs and their lengths
protein_sizes_df = data.frame(
  protein_id = unlist(lapply(prot_info_lst,getName))
  , prot_length = unlist(lapply(prot_info_lst, getLength))) 

protein_sizes_df = protein_sizes_df %>% mutate(protein_id = str_split(protein_id, pattern = "\t", simplify = T)[,1])%>% 
  mutate(protein_id = as.character(protein_id))

rm(prot_info_lst)
gc()

###   2. Count number of proteins per replicon/genome 

prot_number = protein_sizes_df %>% distinct() %>% 
  mutate(replicon =  str_remove(protein_id, pattern="_[0-9]{1,6}$")) %>% 
  group_by(replicon) %>% summarise(protein = n()) %>% ungroup()


###   3. Load and process MMseqs2 default output using mmseqs convertalis qDB tDB rDB result.m8

path_data = "MMseqs2/pp_0523.m8"
names_lst = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen","qstart", "qend" , "sstart", "send" , "evalue", "bitscore")

# Read and calculate alignment coverage from both query and subject sides
mmseqs_out = fread(path_data, col.names = names_lst, header = F) %>%  left_join(protein_sizes_df, by=c("qseqid"="protein_id")) %>% 
  mutate(cov_q = length/prot_length) %>% select(-prot_length) %>% left_join(protein_sizes_df, by=c("sseqid"="protein_id")) %>% 
  mutate(cov_s = length/prot_length) %>% select(-prot_length)

# Parameters for filtering
filtered_df =  mmseqs_out %>% filter(evalue < 1e-4, pident > 0.35 & cov_q > 0.5 & cov_s > 0.5) %>% 
  select(qseqid,sseqid,pident,evalue,bitscore,cov_q,cov_s) %>%
  # set two new columns with name 
  mutate(
    hit_id = paste(qseqid,sseqid, sep = "_")
    ,query = str_remove(qseqid, pattern="_[0-9]{1,6}$")
    ,subject = str_remove(sseqid, pattern="_[0-9]{1,6}$")
    ,pair_id = paste(query, subject, sep = "_")) 

rm(mmseqs_out)
gc()

###   4. Keep bidirectional protein hits

#  Make a table for cross hits
cross_hit_df =  filtered_df %>% select(qseqid, sseqid, subject, query,cross_bitscore=bitscore) %>% 
  mutate(hit_id = paste(sseqid,qseqid, sep = "_")
         ,pair_id = paste(subject,query, sep ="_")) %>% 
  select(hit_id)                                

# Make a semi_join to get only bidirectional hits,           
bidirectional_df = filtered_df %>% semi_join(cross_hit_df, copy = T) 

rm(cross_hit_df, filtered_df)
gc()


# Label the hits according to single hits and multiple hits
bidirectional_df = bidirectional_df %>%
  arrange(qseqid,subject,evalue) %>% group_by(qseqid, subject, pair_id) %>% mutate(n_q=dplyr::n()) %>% 
  ungroup %>% arrange(sseqid,query,evalue) %>% 
  group_by(sseqid, query, pair_id) %>% mutate(n_s = dplyr::n())  %>% ungroup()  

###   5. Separate Single-Hit and Multi-Hit Protein Alignments

# Take single hits
single_hits_df = bidirectional_df %>% filter(n_q ==1 & n_s == 1) %>% group_by(pair_id,query) %>% 
  mutate(single_hit = dplyr::n()) %>% ungroup()

# Handle multiple hits problem

# Multiple hits
multiple_hits_df = bidirectional_df %>% anti_join(single_hits_df, copy = T)

###   6. Resolve multi-hit Cases 

### Make a while loop to get the hits from the multiple hits pool
# Get best hits, remove sequences from the multiple hits table until table is solved => no hits left

# A. Do this one time from the query side 
solved_hits_df_q = data.frame(NULL)
while (length(multiple_hits_df$hit_id)!=0){
  tmp_solv_case = multiple_hits_df %>% group_by(pair_id,subject,qseqid)%>% arrange(desc(bitscore),desc(pident),desc(cov_q))  %>% mutate(n_q=1:n()) %>%
    slice(1) %>% group_by(sseqid,query,pair_id) %>% mutate(n_s=1:n()) %>% arrange(desc(bitscore),desc(pident),desc(cov_s)) %>% 
    slice(1) %>% ungroup()
  solved_hits_df_q = solved_hits_df_q %>% bind_rows(tmp_solv_case)
  multiple_hits_df = multiple_hits_df %>% anti_join(solved_hits_df_q,by=c("qseqid","subject")) %>% 
    anti_join(solved_hits_df_q, by=c("sseqid","query"))   
  tmp_solv_case = data.frame(NULL)
} 

# read in m_hits table
# subject while loop
multiple_hits_df = bidirectional_df %>% anti_join(single_hits_df, copy = T)

# B. Another time from the subject side 
solved_hits_df_s = data.frame(NULL)
while (nrow(multiple_hits_df)!=0){
  tmp_solv_case = multiple_hits_df %>% group_by(sseqid,query,pair_id) %>% arrange(desc(bitscore),desc(pident),desc(cov_q)) %>% 
    mutate(n_s=1:n()) %>% slice(1) %>% group_by(qseqid,subject,pair_id) %>% arrange(desc(bitscore),desc(pident),desc(cov_s)) %>% 
    mutate(n_q=1:n()) %>% slice(1) %>% ungroup()
  solved_hits_df_s = solved_hits_df_s %>% bind_rows(tmp_solv_case)
  multiple_hits_df = multiple_hits_df %>% anti_join(solved_hits_df_s,by=c("sseqid","query")) %>%
    anti_join(solved_hits_df_s, by=c("qseqid","subject"))   
  tmp_solv_case = data.frame(NULL)
} 


###   7. Finalize Best Bidirectional Hits (BBH) table

# Merge single hits and resolved multiple hits, keeping best scoring hit
solved_df= bind_rows(solved_hits_df_q, solved_hits_df_s)%>% select(-n_q, -n_s) %>% 
  group_by(qseqid,subject,pair_id) %>% arrange(desc(bitscore),desc(pident), desc(cov_q)) %>% 
  mutate(n_q=1:n()) %>% group_by(sseqid,query,pair_id) %>% arrange(desc(bitscore),desc(pident)) %>% 
  mutate(n_s=1:n()) %>% filter(n_q ==1 & n_s == 1) %>% arrange(hit_id) %>% ungroup()

# Make the BBH table
BBH_df = bind_rows(single_hits_df,solved_df) 

# Save BBH table to a file
write_tsv(BBH_df, "pp_BBH_df.tsv")


###   8. Calculate weighted Gene Repertoire Relatedness (wGRR) values from both genome sizes

# Include protein counts per genome
score_df =  BBH_df %>%  
  left_join(y = prot_number, by= c("query"="replicon"),copy = T) %>% 
  rename(prot_query = protein) %>%
  left_join(y = prot_number, by= c("subject"="replicon"),copy = T) %>% 
  rename(prot_subject = protein) %>%
  arrange(subject) %>% 
  mutate(ident_s = (pident/prot_subject), 
         ident_q =        (pident/prot_query), 
         edge_id_s = paste(subject,query, sep ="_"),
         edge_id_q = paste(query,subject, sep ="_")) %>%  
  group_by(query,subject) %>% 
  summarize(score_s =  round(sum(ident_s, na.rm = TRUE),3)
            ,score_q =  round(sum(ident_q, na.rm = TRUE),3)
            ,n_similar_proteins =n()
  ) %>% ungroup() 

# Compute wGRR using the genome with fewer proteins
wGRR_df  = score_df  %>% 
  left_join(y = prot_number, by= c("query"="replicon"),copy = T) %>% 
  rename(prot_query = protein) %>%
  left_join(y = prot_number, by= c("subject"="replicon"),copy = T) %>% rename(prot_subject = protein) %>%
  mutate(wGRR = ifelse(prot_subject < prot_query, score_s , score_q) ) %>% 
  select(query, subject, wGRR,n_similar_proteins, prot_subject, prot_query)

# Save wGRR table to a file
write_tsv(wGRR_df, "pp_wGRR_df.tsv")


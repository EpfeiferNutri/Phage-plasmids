# load libraries

library(seqinr)
library(tidyverse)

# # inspect bacterial genome
# 
# ### 204G7 
# 
# # load in genome
# fasta_lst = read.fasta( file = "204G7/204G7_assemblies_flye.fasta"
#                         , as.string = T, forceDNAtolower = F)
# 
# fasta_df = data.frame( contig = unlist(lapply(fasta_lst,getName))
#                        , contig_size = unlist(lapply(fasta_lst,getLength))
#                        , sequence = unlist(lapply(fasta_lst, getSequence, as.string = T)) ) %>%
#   mutate(contig= str_c("nonviral_contig_",contig))
# 
# 
# ### run genomad on the genome, and take a look on the viral contigs
# viral_fasta_lst = read.fasta( file = "204G7/204G7_assemblies_flye_viral.fasta"
#                               , as.string = T, forceDNAtolower = F)
# 
# viral_fasta = data.frame( contig = unlist(lapply(viral_fasta_lst,getName))
#                           , contig_size = unlist(lapply(viral_fasta_lst,getLength))
#                           , sequence = unlist(lapply(viral_fasta_lst,getSequence, as.string = T)) )%>%
#   mutate(contig= str_c("viral_contig_",contig))
# 
# contig_sizes = bind_rows(fasta_df, viral_fasta) %>% select(-sequence)
# 
# 
# ### After mapping, analyze the viral community
# 
# ## profile files, number of reads (total) = 132077
# profile_204G7 = read_tsv("204G7/204G7_rest_long_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#   mutate(contig = str_c("nonviral_contig_", contig))  %>% 
#   bind_rows(read_tsv("204G7/204G7_viral_long_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#               mutate(contig= str_c("viral_contig_",contig))) %>% 
#   left_join(contig_sizes) %>% filter(read_counts > 0) %>%
#   mutate(Host_Viral = ifelse(grepl(contig, pattern="nonviral"),"Host",contig)) %>%
#   group_by(Host_Viral) %>% summarise(contig_size = sum(contig_size), read_counts = sum(read_counts)) %>%
#   # remove viral contigs from sizes
#   mutate(Average_read_cov = round(read_counts/contig_size,5), Strain = "204G7", reads_covered = round(100*sum(read_counts)/132077,1)
#          ,Host_PP = case_when(Host_Viral %in% c("viral_contig_1") ~ "P-P", Host_Viral == "Host" ~ "Host", T ~ "Phage"))
# 
# profile_204G7 = profile_204G7 %>% mutate(Norm_read_cov = Average_read_cov*1/(profile_204G7 %>% filter(Host_Viral == "Host"))$Average_read_cov )
# 
# ### 174J8
# 
# # load in genome
# fasta_lst = read.fasta( file = "174J8/174J8.fasta"
#                         , as.string = T, forceDNAtolower = F)
# 
# fasta_df = data.frame( contig = unlist(lapply(fasta_lst,getName))
#                        , contig_size = unlist(lapply(fasta_lst,getLength))
#                        , sequence = unlist(lapply(fasta_lst, getSequence, as.string = T)) ) %>%
#   mutate(contig= str_c("nonviral_contig_",contig))
# 
# 
# ### run genomad on the genome, and take a look on the viral contigs
# # inspect assembly
# viral_fasta_lst = read.fasta( file = "174J8/174J8_assembly_virus.fna"
#                               , as.string = T, forceDNAtolower = F)
# 
# viral_fasta = data.frame( contig = unlist(lapply(viral_fasta_lst,getName))
#                           , contig_size = unlist(lapply(viral_fasta_lst,getLength))
#                           , sequence = unlist(lapply(viral_fasta_lst,getSequence, as.string = T)) )%>%
#   mutate(contig= str_c("viral_contig_",contig))
# 
# contig_sizes = bind_rows(fasta_df, viral_fasta) %>% select(-sequence)
# 
# 
# ### After mapping, analyze the viral community
# 
# ## profile files, number of reads (total) = 7172960
# profile_174J8 = read_tsv("174J8/174J8_rest_longpac_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#   mutate(contig = str_c("nonviral_contig_", contig))  %>% 
#   bind_rows(read_tsv("174J8/174J8_viral_long_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#               mutate(contig= str_c("viral_contig_",contig))) %>% 
#   left_join(contig_sizes) %>% filter(read_counts > 0) %>%
#   mutate(Host_Viral = ifelse(grepl(contig, pattern="nonviral"),"Host",contig)) %>%
#   group_by(Host_Viral) %>% summarise(contig_size = sum(contig_size), read_counts = sum(read_counts)) %>%
#   # remove viral contigs from sizes
#   mutate(contig_size = ifelse(Host_Viral == "Host", contig_size - sum((contig_sizes %>% filter(!grepl(contig, pattern="nonviral")))$contig_size), contig_size)) %>% 
#   #make average
#   mutate(Average_read_cov = round(read_counts/contig_size,5), Strain = "174J8" , reads_covered = round(100*sum(read_counts)/7172960,1)
#          ,Host_PP = case_when(Host_Viral %in% c("viral_contig_1") ~ "P-P", Host_Viral == "Host" ~ "Host", T ~ "Phage"))
# 
# profile_174J8 = profile_174J8 %>% mutate(Norm_read_cov = Average_read_cov*1/(profile_174J8 %>% filter(Host_Viral == "Host"))$Average_read_cov )
# 
# 
# ### 170D8
# 
# # load in genome
# fasta_lst = read.fasta( file = "170D8/170D8_long_read_assembly.fasta"
#                         , as.string = T, forceDNAtolower = F)
# 
# fasta_df = data.frame( contig = unlist(lapply(fasta_lst,getName))
#                        , contig_size = unlist(lapply(fasta_lst,getLength))
#                        , sequence = unlist(lapply(fasta_lst, getSequence, as.string = T)) ) %>%
#   mutate(contig= str_c("nonviral_contig_",contig))
# 
# 
# ### run genomad on the genome, and take a look on the viral contigs
# # inspect assembly
# viral_fasta_lst = read.fasta( file = "170D8/170D8_long_read_assembly_virus.fna"
#                               , as.string = T, forceDNAtolower = F)
# 
# viral_fasta = data.frame( contig = unlist(lapply(viral_fasta_lst,getName))
#                           , contig_size = unlist(lapply(viral_fasta_lst,getLength))
#                           , sequence = unlist(lapply(viral_fasta_lst,getSequence, as.string = T)) )%>%
#   mutate(contig= str_c("viral_contig_",contig))
# 
# contig_sizes = bind_rows(fasta_df, viral_fasta) %>% select(-sequence)
# 
# 
# ### After mapping, analyze the viral community
# 
# ## profile files, number of reads (total) = 8521077
# profile_170D8 = read_tsv("170D8/170D8_rest_long_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#   mutate(contig = str_c("nonviral_contig_", contig))  %>% 
#   bind_rows(read_tsv("170D8/170D8_viral_long_ab_profile.txt", skip = 11, col_names = c("contig","read_counts")) %>% 
#               mutate(contig= str_c("viral_contig_",contig))) %>% 
#   left_join(contig_sizes) %>% filter(read_counts > 0) %>%
#   mutate(Host_Viral = ifelse(grepl(contig, pattern="nonviral"),"Host",contig)) %>%
#   group_by(Host_Viral) %>% summarise(contig_size = sum(contig_size), read_counts = sum(read_counts)) %>%
#   # remove viral contigs from sizes
#   mutate(contig_size = ifelse(Host_Viral == "Host", contig_size - sum((contig_sizes %>% filter(!grepl(contig, pattern="nonviral")))$contig_size), contig_size)) %>% 
#   #make average
#   mutate(Average_read_cov = round(read_counts/contig_size,5), Strain = "170D8", reads_covered = round(100*sum(read_counts)/8521077,1)
#          ,Host_PP = case_when(Host_Viral %in% c("viral_contig_16") ~ "P-P", Host_Viral == "Host" ~ "Host", T ~ "Phage"))
# 
# profile_170D8 = profile_170D8 %>% mutate(Norm_read_cov = Average_read_cov*1/(profile_170D8 %>% filter(Host_Viral == "Host"))$Average_read_cov )
# 
# 
# ## combine and plot
# profiles_plot_df = bind_rows(profile_170D8, profile_174J8, profile_204G7) 
# 
# profiles_plot_df %>% write_tsv("read_count_table.tsv")

read_count_plot_tbl = read_tsv("read_count_table.tsv")

read_count_plot_tbl %>%
  filter(Host_PP != "Host") %>%
  ggplot(aes(x=Norm_read_cov, y=Host_Viral, fill = Host_PP)) + 
  geom_col(width = 0.6, col = "black") + theme_bw()+
  scale_fill_manual(values = c("orange","royalblue"))+scale_x_log10()+
  ylab("")+xlab("Normalized read coverage (au)")+
  theme(text = element_text(size = 14), strip.background = element_blank(), legend.title = element_blank(), legend.position = "top")+
  facet_wrap(~Strain,nrow = 3, scales = "free" )




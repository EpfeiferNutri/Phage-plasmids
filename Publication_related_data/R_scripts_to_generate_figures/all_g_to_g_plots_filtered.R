
########## Script for generating genome-to-genome comparative plots for selected genomes

# required libraries

library(data.table) 
library(gggenomes)  
library(ggplot2)   
library(dplyr)     
library(ggnewscale)
library(tidyverse) 

############ Reading the data tables required for generating the figures

# Genome sequences used for plotting 
genomes = fread("Publication_related_data/Supplementary_data/g2g_plot_tables/genomes_used_gggenomes.tsv")

# Protein information necessary for gggenomes plots (coordinates, protein IDs, etc.)
protein_info = fread("Publication_related_data/Supplementary_data/g2g_plot_tables/protein_info_used_gggenomes.tsv")

# Best Bi-directional Hit (BBH) protein-protein links
# Note: seq_id == query, seq_id2 == subject
BBH_links = fread("Publication_related_data/Supplementary_data/g2g_plot_tables/BBH_links_used_gggenomes.tsv")

# Table with wGRR (weighted Gene Repertoire Relatedness) values for selected genomes 
wGRR_info = fread("/Publication_related_data/Supplementary_data/g2g_plot_tables/wGRR_info_used_gggenomes.tsv")

############ Defining the genome-to-genome comparative plotting function

plot_gggenomes = function(genomes_list_, protein_info_, links_) {
  
  gggenomes(seqs = genomes_list_, genes = protein_info_, links = links_) +
    geom_seq() +        
    geom_gene(aes(fill = is_persistent)) +   # color genes by persistence status
    scale_fill_manual(values = c("limegreen", "grey80")) +
    
    new_scale_fill() +  
    
    geom_link(aes(fill = pident), color = alpha("white", alpha = 0), offset = 0.08) +
    scale_fill_gradient2(
      limits = c(0.35, 1), 
      breaks = c(0.35, 0.5, 0.65, 0.8, 0.95),
      low = alpha("white", alpha = 0.0),
      midpoint = 0.5, 
      mid = alpha("tomato", alpha = 1), 
      high = alpha("darkgreen", alpha = 1)) +
    
    guides(colour = guide_colorbar(), direction = "horizontal") +   # colorbar for homology identity
    geom_bin_label(size = 8, expand_left = 0.2) +
    
    theme(
      axis.text.x = element_text(size = 20),
      legend.title = element_blank(),
      # legend.box = "horizontal",
      legend.position = "right") +
    
    scale_x_bp(suffix = "b", sep = " ")
  
}

############# Defining the processing function: prepare data for a specific Phage-Plasmid group

process_pp_type = function(PP_type, genome_list_i) {
  
  # Subset genomes and order according to user-defined list
  genomes_i = genomes %>% filter(seq_id %in% genome_list_i)%>%
    mutate(seq_id = factor(seq_id, levels = genome_list_i)) %>% arrange(seq_id)
  
  # Subset BBH links and protein information for selected genomes
  links_i = BBH_links %>% filter(seq_id %in% genomes_i$seq_id)
  protein_info_i = protein_info %>% filter(seq_id %in% genomes_i$seq_id)
  
  # Generate gggenomes plot
  plot_gggenomes(genomes_list_ = genomes_i, protein_info_ = protein_info_i, links_ = links_i)

}

############ Generating the plots for each figure and dataset

################# Examples of hits (low confidence) with atypical sizes: 
# --------------- pMT1 long genomes ----------------
genome_list_pMT1 = c("NZ_CP045146", "NZ_CP010248")
process_pp_type(PP_type = "pMT1", genome_list_i = genome_list_pMT1)
# --------------- pKpn validation ----------------
genome_list_pKpn = c("NZ_CP017778",'NZ_CP028550', 'NZ_CP027149')
process_pp_type(PP_type = "pKpn", genome_list_i = genome_list_pKpn)

################# cp32 elements with high, medium and low confidence levels:
genome_list_cp32 = c("NZ_CP115597", "NZ_CP114697", "NZ_CP073149", "NZ_CP019758", "NZ_CP019926")
process_pp_type(PP_type = "cp32", genome_list_i = genome_list_cp32)

################# Detection of P1-like plasmids (false-positives) by tyPPing:
genome_list_P1_1_1 = c('NC_005856', 'NZ_CP010146', "NZ_LR883966", "NZ_CP061756")
process_pp_type(PP_type = "P1_1", genome_list_i = genome_list_P1_1_1)

################# N15-like P-Ps with atypical sizes:
genome_list_N15 = c("NZ_CP064249", "NZ_CP070567", "NZ_CP082016")
process_pp_type(PP_type = "N15", genome_list_i = genome_list_N15)

################# P1_1 low confidence with too-long size:
genome_list_P1_1_2 = c("NZ_CP081715",'NZ_CP103719', 'NC_005856')
process_pp_type(PP_type = "P1_1", genome_list_i = genome_list_P1_1_2)

################# P1_2 case (false positive) predicted by tyPPing in 05/23:
genome_list_P1_2 = c("NZ_CP100729", "NZ_CP032492")
process_pp_type(PP_type = "P1_2", genome_list_i = genome_list_P1_2)

################# Cases (low confidence) related to AB_1, detected by MM-GRC and not by tyPPing:
genome_list_AB = c( 'NZ_CP107600',                            # virus, +tyPPing
                    "NZ_CP069498",                            # provirus, -tyPPing
                    "NZ_CP102935", "NZ_CP102763")             # provirus, -tyPPing, -vcontact
process_pp_type(PP_type = "AB", genome_list_i = genome_list_AB)

################# SSU5_pHCM2-like P-P not detected by tyPPing :
genome_list_SSU5 = c("NZ_CP110081", "NZ_CP077310")
process_pp_type(PP_type = "SSU5_pHCM2", genome_list_i = genome_list_SSU5)

################# pSLy3-like P-P not detected by tyPPing:
genome_list_pSLy3 = c("NZ_CP031294", "NZ_CP074043")
process_pp_type(PP_type = "pSLy3", genome_list_i = genome_list_pSLy3)


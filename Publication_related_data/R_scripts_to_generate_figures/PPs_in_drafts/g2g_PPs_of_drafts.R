# libraries
library(tidyverse)
library(gggenomes)
library(readxl)

#### make a plot of the genomes using gggenomes

   # load in genome table (sizes, accessions, seq id, etc.)
    pps_df = read_tsv("data_PPs_induction/g2g_plots/pps_table.tsv")
 
   # PP gene tbl
   genes_table = read_tsv("data_PPs_induction/g2g_plots/PP_genes_table.tsv")
 
 ### BBH tbl
   BBH_tbl = read_tsv("data_PPs_induction/g2g_plots/BBH_tbl.tsv")
 
  # wGRR, not needed for plot, just for inspection
   wGRR = read_tsv("data_PPs_induction/g2g_plots/wGRR_tbl.tsv")
    
   # select genomes to plot
   plot_df = pps_df %>% 
    #filter(seq_id %in% c("174J8_P1_2","NC_005856")) %>%
    mutate(seq_id=NCBI_accession) %>%
    select(seq_id,length=size,replicon,NCBI_accession) %>%
    mutate(order = case_when(seq_id == "174J8_P1_2" ~ 1
            ,seq_id == "NC_005856" ~ 2,seq_id == "170D8_P1_1" ~ 3
            ,seq_id == "NC_050154" ~ 4,seq_id == "204G7_N15" ~ 5, seq_id == "NC_001901" ~ 6)) %>% arrange(order)
   
   
## make the link table
 links_df = BBH_tbl %>% 
   left_join(genes_table %>% select(qseqid=gene, start, end, strand1=strand, seq_id ), by = "qseqid") %>%
   left_join(genes_table %>% select(sseqid=gene, start2=start, end2=end,strand2=strand, seq_id2=seq_id), by = "sseqid") %>%
  mutate(length=end-start, strand = ifelse(strand1==strand2, "+","-")) %>% filter(qseqid!=sseqid,start<end & start2 < end2) %>% 
  select(feat_id = qseqid, feat_id2=sseqid, evalue,pident, bitscore, length ,
         start, end, strand, seq_id, start2, end2,seq_id2)
   
  ### plot make a simple plot
  gggenomes(genes = genes_table, seqs = plot_df) %>%
    add_feats() %>%
    add_links(links_df) %>%
    flip_seqs() +
    geom_seq( col = "grey50", size = 1) +
  geom_link(aes(fill = pident), color = alpha("white", alpha =0), offset = 0.05, position ="identity") +
  scale_fill_gradient2(limits = c(0.35,1), breaks = c(0.35, 0.5, 0.65,0.8,0.95)
                       , low = alpha("white", alpha = 0.0), midpoint = 0.5, mid = alpha("tomato", alpha =1) , high = alpha("steelblue", alpha = 1))+
  guides(colour=guide_colourbar(),direction = "horizontal")+
  geom_gene(fill = "grey" ,  size = 3 , position = "identity", show.legend = T) +
  geom_bin_label(size = 4, expand_left = 0.2) +
  theme(axis.text.x=element_text(size=12)
        , legend.title = element_blank(), legend.box = "horizontal"
  , legend.position = "right") +
  scale_x_bp(suffix = "b", sep = " ") 
  

  
  
  
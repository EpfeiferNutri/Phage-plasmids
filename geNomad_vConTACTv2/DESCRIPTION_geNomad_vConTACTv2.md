# Introduction

**geNomad** is a tool for classifying and annotating mobile genetic elements as phages, integrated phages or plasmids. It uses gene content analysis and deep learning approaches. It is not specifically designed for P-P detection, but it can serve as a useful screen. We tested geNomad by looking for phage sequences in a plasmids database.

**vConTACT v2** is a viral classification tool that groups seqeuneces with a reference dataset using gene-sharing networks. We detected P-Ps by geNomad, and used vConTACT v2 to cluster them with a reference dataset of typed P-Ps. vConTACT v2 effectively grouped these sequences with a comparabale accuracy as our previous approach (MM-GRC, and tyPPing).

**Quick guide:**

1.  Run geNomad on a set of sequences that are described to be plasmids.

2.  Interpretation of geNomad output:

    -   If classified as phage → putative P-P

    -   If classified as prophage → reject as P-Ps are not integrated prophages, or undergo further evaluation using complementary tools

    -   If classified as plasmid → reject, since it is likely a true plasmid (not a P-P)

3.  Run vConTACT v2 on the set of putative P-Ps from the previous step combined with a reference P-P dataset. You can use the 1416 P-Ps from PMID: 38378896, or (recommended) filter this datatst for P-Ps with high confidence (assigned by tyPPing).

4.  Interpret clustering results:

    -   If a candidate P-P clusters with a typed P-P → P-P of the corresponding type.

    -   If it clusters elsewhere or remains unclustered → unknown or not a P-P

    -   If it clusters with multiple typed P-Ps → likely a P-P, a co-integration or recombination product

# Usage

### **Tool setup & Installation**

For installation instructions and detailed user guides, see [geNomad documentation](https://portal.nersc.gov/genomad/index.html) and [vConTACTv2 documentation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/).

### Input data **requirements**

**geNomad:**

-   Nucleotide FASTA file 

**vConTACT v2:**

-   A protein sequences FASTA file 
-   A "protein-to-genome" mapping file (`.tsv` or `.csv` format).
-   **For P-P detection**: include a reference P-P dataset in the protein and protein-to-genome tables.
  Use --db flag "None"

### **Running**

**Commands** to execute geNomad and vConTACT v2:

```{bash, eval = FALSE}

genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE

vcontact2 --raw-proteins [proteins file] --proteins-fp [gene-to-genome mapping file] --db None --output-dir [target output directory]
```

# Interpretation and general recommendations

-   geNomad is a fast and scalable tool for screening large plasmid datasets and can effectively detect P-Ps. However, it may produce conflicting predictions (P-Ps detected as integrated prophages or plasmids) and for P-Ps that may lack clear phage signatures (such as cp32). For unsure candidates, cross results of geNomad with predictions of tyPPing and MM-GRC. 

-   vConTACT v2 accuracy depends on the quality of the used reference dataset. To improve performance, combine vConTACT v2 with additional tools like MM-GRC and tyPPing.

# Results summary on 05/23 dataset

```{bash eval=FALSE, include=FALSE}

geNomad_vConTACTv2/
├── DESCRIPTION_geNomad_vConTACTv2.md  
├── geNomad/
│   ├── 0523_plasmids_for_geNomad.xlsx                 # Input list of plasmid sequences
│   ├── geNomad_plasmid_summary_on_0523_plasmids.tsv   # Predicted plasmids (raw geNomad output)
│   ├── geNomad_virus_summary_on_0523_plasmids.tsv     # Predicted phages and prophages (raw geNomad output)
│
├── vConTACTv2/
│   ├── 0523_plasmids_ref_PPs_for_vConTACT.csv         # Input genomes for vConTACT v2 (phages by geNomad (target) + reference P-Ps)
│   ├── reference_0321_PP_list_for_vConTACT2.xlsx      # Reference P-P dataset with assigned P-P types
│   ├── genome_by_genome_overview.csv                  # Pairwise genome clustering overview (raw vConTACT output)
│   ├── c1.ntw                                         # Gene-sharing network file (raw vConTACT output)
├── vConTACT2_PP_list_in_0523_plasmids.xlsx            # Final output: predicted P-Ps with assigned types
```

**Step 1: running geNomad**

-   As input dataset we processed 38,051 plasmid sequences from non-redundant NCBI RefSeq, collected in May 2023 (referred to as "05/23").

-   geNomad (end-to-end) was executed in default mode with the plasmid sequences listed in `0523_plasmids_for_geNomad.xlsx`. It used the geNomad database v1.5 and was executed on 25 CPUs on the Migale computational cluster [https://migale.inrae.fr/], with a total runtime of 6 hours and 17 minutes.

-   Standard geNomad output summary tables (from `<prefix>_summary/` folder) include replicons from 05/23 dataset predicted as plasmids (`geNomad_plasmid_summary_on_0523_plasmids.tsv`) and ones predicted as phages or prophages (`geNomad_virus_summary_on_0523_plasmids.tsv)`.

-   Replicons listed in `geNomad_virus_summary_on_0523_plasmids.tsv` and classified as not prophages (topology != "Provirus") were selected as potential P-Ps.

-   Detection summary for 05/23 dataset (see Table S1):

    -   The virus summary table included 2,333 viral predictions, with 1,714 classified as phages (putative P-Ps) and 619 as prophages.

    -   Among these, 804 new P-Ps are novel in 05/23, (excluding 910 already known in 03/21), and 407/804 overlap with the well-defined P-P types.

    -   geNomad recovered between 69.2% and 79.5% of known P-Ps, i.e. 20-30% were misclassified as plasmids (such as cp32) or as prophages (many of AB_1).

**Step 2: running vConTACT v2**

We used vConTACT v2 to cluster putative P-Ps and assign types with a reference dataset (1416 P-P from PMID: 38378896).

-   As input dataset we processed 1,714 putative P-Ps (`from 0523_plasmids_ref_PPs_for_vConTACT.csv`).

-  vConTACT v2 required 13 minutes (default settings, 25 CPUs on the Migale cluster)

-   Output files from vConTACT v2 include clustering information for each genome (`genome_by_genome_overview.csv`) and gene-sharing network results (`c1.ntw`).

-   Reference P-P list  with 1416 typed P-Ps `reference_0321_PP_list_for_vConTACT2.xlsx`.

-   R code used for assigning types to putative P-Ps based on clustering with reference P-Ps:

```{r}
# P-P reference list + geNomad virus-positive replicons in 05/23

# Load metadata and filter for new replicons (present in 05/23 but absent in 03/21)
new_replicons = read_excel("Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx") %>% 
  filter(`in 05/23` == "Yes", `in 03/21` == "No") %>%
  select(`NCBI accession`)

reference_0321_PP_list_for_vConTACT2 = read_excel("geNomad_vConTACTv2/vConTACTv2/reference_0321_PP_list_for_vConTACT2.xlsx")

genome_by_genome_overview_0523 = fread('geNomad_vConTACTv2/vConTACTv2/genome_by_genome_overview.csv')

# Join clustering results with reference P-P types
PPs_only_genomad_Summary_0523_vcontact = genome_by_genome_overview_0523 %>% 
  left_join(reference_0321_PP_list_for_vConTACT2, by=c("Genome"="NCBI accession")) %>%
  select(Genome, `VC Status`, VC, `VC Subcluster`,`P-P type`)

# Process clustering information to assign P-P types
PPs_only_genomad_Summary_0523_vcontact_processed = PPs_only_genomad_Summary_0523_vcontact %>% 
  # `Virus cluster` is the most general cluster 
  mutate(`Virus cluster` = case_when(`VC Status` == "Singleton" ~ str_c(row_number()," Singleton"), 
                                     `VC Status` == "Outlier" ~ str_c(row_number()," Outlier"),
                                     `VC Status` == "Clustered" ~ str_remove(`VC`, pattern="_[0-9]+$"), 
                                     str_detect(`VC Status`, "Overlap") ~ str_remove_all(`VC Status`, pattern = "^Overlap \\(|\\)$"), 
                                     `VC Status` == "Clustered/Singleton" ~ str_remove(`VC`, pattern="_[0-9]+$"),
                                     T ~ NA),
         `P-P type vConTACT2` = case_when(!is.na(`P-P type`) ~ `P-P type`, 
                                          T ~ NA)) %>%
  # prediction by vcontact;  if conflicting types in one cluster --> all types kept
  group_by(`Virus cluster`) %>%
  mutate(`P-P type vConTACT2` = if (all(is.na(`P-P type`))) NA else paste(unique(na.omit(`P-P type`)), collapse = "; ")) %>%
  ungroup() %>%
  filter(!is.na(`P-P type vConTACT2`)) 

# Extract only newly identified putative P-Ps
putative_new_PPs = PPs_only_genomad_Summary_0523_vcontact_processed %>%
  semi_join(new_replicons, by=c("Genome"="NCBI accession")) %>%
  select(Genome, `P-P type vConTACT2`)
```

-   `vConTACT2_PP_list_in_0523_plasmids.xlsx`: Predicted and typed P-Ps of the 05/23 dataset (detected and typed with geNomad and vConTACT2).

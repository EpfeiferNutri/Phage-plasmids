# Introduction

**geNomad** [REF] is a tool for classifying and annotating mobile genetic elements as phages or plasmids by combining gene content analysis with deep learning. Although it is not specifically designed for **phage-plasmid (P-P)** detection, it can serve as a useful first screen for potential P-Ps by identifying plasmid sequences with phage-like features.

**vConTACT v2** [REF] is a viral classification tool that clusters sequences based on gene-sharing networks. When plasmids identified as phage-like by geNomad (putative P-Ps) are analyzed together with a reference dataset of known P-Ps, vConTACT v2 can effectively group them into the known types.

This document describes the procedure for using geNomad and vConTACT v2 on plasmid sequences and explains how to interpret the results in the context of P-P detection and typing.

**Quick interpretation guide:**

1.  Run geNomad on a set of plasmid sequences.

2.  Interpret the geNomad output based on the assigned classification:

    -   If classified as phage → candidate phage-plasmid (P-Ps)

    -   If classified as prophage → further evaluation recommended using complementary tools (should be excluded unless supported by additional evidence)

    -   If classified as plasmid → likely a true plasmid (not a P-P)

3.  Run vConTACT v2 on the set of putative P-Ps from the previous step combined with a reference P-P dataset.

4.  Interpret clustering results:

    -   If a candidate P-P clusters with known P-Ps → consider it a putative P-P of the corresponding type.

    -   If it clusters elsewhere or remains unclustered → likely not a P-P, or it may be too diverse from the known types.

    -   If multiple reference P-Ps from different types cluster together → the assigned type is identifined as a concatenation of all relevant types.

# Usage

### **Tool setup & Installation**

For installation instructions and detailed user guides, see [geNomad documentation](https://portal.nersc.gov/genomad/index.html) and [vConTACTv2 documentation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/).

### Input data **requirements**

**geNomad:**

-   Any FASTA file containing nucleotide sequences.

**vConTACT v2:**

-   A FASTA file of amino acid sequences.
-   A "gene-to-genome" mapping file (`.tsv` or `.csv` format).
-   **For P-P detection**: a reference dataset of known P-Ps.

### **Running**

**General commands** to execute geNomad [REF] and vConTACT v2 [REF]:

```{bash, eval = FALSE}

genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE

vcontact2 --raw-proteins [proteins file] --proteins-fp [gene-to-genome mapping file] --output-dir [target output directory]
```

# Interpretation and general recommendations

-   geNomad is a fast and scalable tool for screening large plasmid datasets and can effectively detect P-Ps with strong phage signatures. However, it may produce conflicting predictions, especially for integrated or hybrid elements, and for P-Ps that may lack clear phage signatures (such as cp32).

-   For reliable results, combine geNomad with MM-GRC to detect P-P singletons and P-Ps from diverse communities.

-   vConTACT v2 is not tailored for P-P detection, but it can be used for type assignment when applied to target and reference datasets combined. Its accuracy depends on the quality of both datasets.

-   For best performance, combine vConTACT v2 with additional tools like MM-GRC to improve type assignment and reduce false results.

# Results summary on 05/23 dataset

**Folder structure:**

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

We screened the large plasmid dataset to detect putative P-Ps with strong phage signatures.

-   As input dataset we processed 38,051 plasmid sequences from non-redundant NCBI RefSeq, collected in May 2023 (referred to as "05/23").

-   geNomad (end-to-end) was executed in default mode with the plasmid sequences listed in `0523_plasmids_for_geNomad.xlsx`. It used the geNomad database v1.5 and was executed on 25 CPUs on the Migale computational cluster [REF], with a total runtime of 6 hours and 17 minutes.

-   Standard geNomad output summary tables (from `<prefix>_summary/` folder) include replicons from 05/23 dataset predicted as plasmids (`geNomad_plasmid_summary_on_0523_plasmids.tsv`) and ones predicted as phages / prophages (`geNomad_virus_summary_on_0523_plasmids.tsv)`.

-   Replicons listed in `geNomad_virus_summary_on_0523_plasmids.tsv` and classified as phages (topology != "Provirus") were selected as potential P-Ps.

-   Detection summary for 05/23 dataset (see Table S1):

    -   The virus summary table included 2,333 viral predictions, with 1,714 classified as phages (putative P-Ps) and 619 as prophage-containing replicons.

    -   Among these, 804 new putative P-Ps were identified (after excluding 910 known P-Ps), and 407/804 were assigned to well-defined P-P types.

    -   geNomad recovered between 69.2% and 79.5% of known P-Ps, meaning some were misclassified as plasmids (such as cp32) or as prophages (such as AB_1).

**Step 2: running vConTACT v2**

In the next step, we used vConTACT v2 to cluster putative P-Ps and assign types based on gene-sharing with a reference dataset (1416 P-P from [REF]).

-   As input dataset we processed 2,267 sequences (`from 0523_plasmids_ref_PPs_for_vConTACT.csv`), including all putative P-Ps classified as phages by geNomad and the reference P-P dataset.

-   Run details: vConTACT v2 ran with default settings on 25 CPUs on the Migale computational cluster with a total runtime of 13 minutes.

-   The raw output files from vConTACT v2 include clustering information for each genome (`genome_by_genome_overview.csv`) and gene-sharing network results (`c1.ntw`).

-   Reference 03/21 P-P list which was used for vConTACT v2 type assignment is available in `reference_0321_PP_list_for_vConTACT2.xlsx`.

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

-   The final list of predicted P-Ps from the 05/23 dataset (based on combined geNomad and vConTACT2 results) is available in `vConTACT2_PP_list_in_0523_plasmids.xlsx`.

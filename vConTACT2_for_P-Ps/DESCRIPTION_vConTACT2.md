# Introduction

**vConTACT v2** [REF] is a viral classification tool that organizes viral sequences into clusters based on gene-sharing networks.

While it is designed for viral taxonomy, it can be effectively applied to identify and type phage-plasmids (P-Ps) when used in combination with a curated reference dataset.

This document describes the procedure for using vConTACT v2 on plasmid sequences and explains how to interpret the results in the context of P-P detection and typing.

**Quick interpretation guide:**

0.  Run geNomad or RF classification models from MM-GRC on a set of plasmid sequences (for details, see the corresponding documentation).

1.  Run vConTACT v2 on the **set of putative P-Ps combined with reference P-P dataset** (provided in this study).

2.  Interpret clustering results:

    -   If a plasmid clusters with known P-Ps → consider it a putative P-P of the corresponding type.

    -   If it clusters elsewhere or remains unclustered → likely not a P-P.

# Usage

## **Tool Setup & Installation**

For installation instructions and detailed setup, see [vConTACTv2 documentation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/).

## Input data **requirements**

-   A FASTA file of amino acid sequences.
-   A "gene-to-genome" mapping file (`.tsv` or `.csv` format).
-   **For P-P detection**: a curated reference dataset of known P-Ps.

## **Running** vConTACT v2

**General command** to execute vConTACT v2 [REF]:

```{bash, eval = FALSE}
vcontact2 --raw-proteins [proteins file] --proteins-fp [gene-to-genome mapping file] --output-dir [target output directory]
```

**Example run**

-   **Dataset used in this study [REF]:**\
    All sequences were retrieved from the non-redundant NCBI RefSeq database [REF]. The dataset (referred to as *"05/23"*) contains **38,051 plasmid sequences**, collected in **May 2023**.
-   For this analysis, we selected:
    -   All **putative P-Ps** predicted as phages by **geNomad** (see *geNomad_for_P-Ps*).
    -   A **reference set** of **1,416 known P-Ps** from [REF].
-   **vConTACT** was executed in default mode using the 2,267 sequences listed in `0523_plasmids_ref_PPs_for_vConTACT.csv`.

```{r message=FALSE, warning=FALSE, paged.print=FALSE}
library(readxl)
library(data.table)
library(tidyverse)

plasmids_0523_ref_PPs_for_vConTACT = fread("vConTACT2_for_P-Ps/0523_plasmids_ref_PPs_for_vConTACT.csv")
plasmids_0523_ref_PPs_for_vConTACT
```

**Run details:**

-   **Hardware**: 25 CPUs
-   **Total runtime**: 13 minutes

# Results

## Standard vConTACTv2 output summary

The primary output files from vConTACT v2 include clustering and gene-sharing network results.

-   **`genome_by_genome_overview.csv`**: Contains clustering information for each genome, confidence scores and related metrics.

**Genome by genome overview for 05/23 dataset**:

```{r}
genome_by_genome_overview_0523 = fread('vConTACT2_for_P-Ps/genome_by_genome_overview.csv')

genome_by_genome_overview_0523
```

-   **`c1.ntw`**: Represents the gene-sharing network.

**Network data for 05/23 dataset**:

```{r message=FALSE, warning=FALSE}
significance_scores_0523 = fread("vConTACT2_for_P-Ps/c1.ntw") %>% 
  rename(Genome_1 = V1, Genome_2 = V2, Significance_score = V3) 

significance_scores_0523
```

## Classification strategy for putative P-Ps

-   Sequences from a plasmid dataset detected as phages by geNomad should be clustered together with the reference set of known P-Ps [REF].

-   Any plasmid that clustered with a known P-P is considered a putative P-P and assigned the same type.

-   If multiple reference P-Ps from different types cluster together, the assigned type identifined as a concatenation of all associated types.

### vConTACTv2 summary processing for 05/23 dataset

-   Reference 03/21 P-P list used for vConTACT v2 (contains 1,416 known P-Ps) [REF].

```{r}
reference_0321_PP_list_for_vConTACT2 = read_excel("vConTACT2_for_P-Ps/reference_0321_PP_list_for_vConTACT2.xlsx")

reference_0321_PP_list_for_vConTACT2
```

-   Assigning types to putative P-Ps based on clustering with reference P-Ps.

This step integrates geNomad positive replicons from 05/23 with the vConTACT clustering results and assigns P-P types by matching to the reference dataset.

```{r}
# P-P reference list + geNomad virus-positive replicons in 05/23

# Load metadata and filter for new replicons (present in 05/23 but absent in 03/21)
new_replicons = read_excel("Study_associated_data/Supplementary_data/S1_All_MGE_information_table.xlsx") %>% 
  filter(`in 05/23` == "Yes", `in 03/21` == "No") %>%
  select(`NCBI accession`)

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

-   P-P list for 05/23 vConTACT2 predictions

```{r}
pps_0523_by_vcontact = read_excel("vConTACT2_for_P-Ps/vConTACT2_PP_list_in_0523_plasmids.xlsx")

pps_0523_by_vcontact
```

## Detection summary for 05/23 plasmids

-   Total number of typed putative P-Ps: 614

-   Number of new P-Ps clustered into well-defined types: xxxxx

-   Consistently typed P-Ps (detected by vConTACT v2, tyPPing, and MM-GRC): xxxxx

-   Inconsistent cases (conflicting classifications across tools): xxxxx

## Interpretation and general recommendations

-   While vConTACT v2 is not specifically designed to detect and cluster P-Ps, it performs well when used with curated datasets.

-   Accurate type assignment depends heavily on high-quality reference and target datasets

-   In this analysis, we first applied a geNomad to preselect candidate P-Ps for clustering.

-   For best results, combine vConTACT v2 with other approaches such as geNomad and MM-GRC to improve both recall and precision.

-   Be aware that sequence length can impact results (may be misclassified or unclustered), especially for unusually short or long sequences.

# References

-   [tyPPing]
-   [vConTACTv2 documentation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/)

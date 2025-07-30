# Introduction

**geNomad** is a computational tool for classification and annotation of mobile genetic elements. It integrates gene content analysis with deep neural network models to identify sequences corresponding to plasmids and viruses.

Although geNomad is not specifically tailored for detecting **phage-plasmids (P-Ps)**, its phage classification capabilities can be applied to plasmid datasets to identify potential P-P candidates.

This document outlines the workflow for applying geNomad to a set of plasmid sequences and provides guidance on how to interpret the results in the context of P-P detection.

**Quick interpretation guide:**

1.  Run geNomad on a set of plasmid sequences.

2.  Interpret the output based on the assigned classification:

    -   If classified as phage → candidate phage-plasmid

    -   If classified as prophage → further evaluation recommended using complementary tools (a potential recombination event)

    -   If classified as plasmid → likely a true plasmid (not a P-P)

# Usage

## **Tool setup & Installation**

For installation instructions and detailed setup, see [geNomad documentation](https://portal.nersc.gov/genomad/index.html).

## Input data **requirements**

-   **Accepted input:** any FASTA file containing nucleotide sequences (e.g., isolate genomes, metagenomes, or metatranscriptomes) [REF].

## **Running** geNomad

**General command** to execute geNomad [REF]:

```{bash, eval = FALSE}
genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
```

**Example run**

-   **Dataset used in this study [REF]:**\
    All sequences were retrieved from the non-redundant NCBI RefSeq database [REF]. The dataset (referred to as *"05/23"*) contains **38,051 plasmid sequences**, collected in **May 2023**.
-   **geNomad (end-to-end)** was executed in default mode with the plasmid sequences listed in `0523_plasmids_for_geNomad.xlsx`.

```{r message=FALSE, warning=FALSE}
library(readxl)
library(data.table)
library(tidyverse)

plasmids_0523_for_geNomad = read_excel("geNomad_for_P-Ps/0523_plasmids_for_geNomad.xlsx")
plasmids_0523_for_geNomad 
```

**Run details:**

-   **Version**: 1.7.4
-   **Database**: geNomad DB v1.5
-   **Hardware**: 25 CPUs
-   **Total runtime**: 6 hours 17 minutes

# Results

## Standard geNomad output summary tables

The *\<prefix\>\_summary* directory contains files that summarize the results that were generated across the pipeline.

-   Replicons from 05/23 dataset predicted as plasmids:

```{r}
geNomad_plasmid_summary_on_0523_plasmids = fread("geNomad_for_P-Ps/geNomad_plasmid_summary_on_0523_plasmids.tsv")

geNomad_plasmid_summary_on_0523_plasmids
```

-   Replicons from 05/23 dataset predicted as phages / prophages:

```{r}
geNomad_virus_summary_on_0523_plasmids = fread("geNomad_for_P-Ps/geNomad_virus_summary_on_0523_plasmids.tsv")

geNomad_virus_summary_on_0523_plasmids
```

## Classification strategy for detecting putative P-Ps in 05/23 dataset :

We applied the following criteria to detect putative P-Ps from the geNomad results:

-   Replicons classified as phages were considered putative P-Ps.

-   Replicons classified as integrated prophages were excluded unless supported by additional evidence.

-   Replicons predicted as plasmids were retained as such unless flagged by other tools.

### Processing the geNomad's output for 05/23 dataset

-   Identifying potential P-Ps (extract replicons classified as phages, excluding those labeled as prophages)

```{r}
phages_by_geNomad = geNomad_virus_summary_on_0523_plasmids %>%  
  select(-fdr, -genetic_code, -coordinates) %>% 
  filter(topology != "Provirus")

phages_by_geNomad
```

-   Identifying and summarizing prophage-carrying replicons

```{r}
prophages_by_geNomad = geNomad_virus_summary_on_0523_plasmids %>% filter(topology == "Provirus")

prophage_containing_replicons_by_geNomad = prophages_by_geNomad %>% mutate(replicon = str_remove(seq_name,  "\\|provirus.*$")) %>% 
  group_by(`replicon`) %>% summarize(count_prophages = n()) 

prophage_containing_replicons_by_geNomad
```

## Detection summary for 05/23 dataset

-   **Total viral hits in geNomad output**: 2333 sequences

    -   1,714 predicted as **phages**

    -   619 predicted as **prophage-containing**

-   **Putative new P-Ps identified**: 804 sequences (derived by subtracting 910 previously known P-Ps (from 03/21 dataset) from the 1,714 phage-classified sequences)

    -   **407/804** belong to well-defined P-P types (Table S1).

-   geNomad recovered **69.2% to 79%** of previously predicted P-Ps.

## Interpretation and general recommendations

-   **geNomad** is a fast and scalable tool for initial screening of P-Ps in large **plasmid** datasets.

-   It is effective for identifying P-Ps **with strong phage signatures**, but results should be interpreted with caution.

-   geNomad frequently causes conflicting P-P predections, if solely used (especially for integrated or hybrid elements).

-   Use geNomad **on confident plasmid sequences** in combination with MM-GRC to detect P-P singletons and P-Ps from diverse communities.

-   Double-check sequences that are borderline or with inconsistent classifications across tools (such as MM-GRC or vConTACT v2).

# References

-   [tyPPing]
-   [geNomad documentation](https://portal.nersc.gov/genomad/index.html)

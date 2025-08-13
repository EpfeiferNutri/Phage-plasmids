# Introduction

Multi-model gene repertoire clustering, or short **MM-GRC,** was designed to detect and type P-Ps in phage and plasmid datasets. It detects P-Ps of well-defined and other types (including singletons).

**MM-GRC** needs to be run on plasmid and phage sequences separetly. 

**MM-GRC workflow**

<img width="2074" height="794" alt="image" src="https://github.com/user-attachments/assets/a7a4ff2d-ceb4-440f-a184-0fe5afec85c3" />

1.  Annotate target proteins of phages and plasmids using specific HMMs (phage sequences with plasmids profiles and vice versa for plasmid sequences).

2.  
    A) Use random forest classifiers to detect P-Ps in plasmid databases. Consider only elements with genome sizes between 10 kb and 300 kb.
    B) Consider phage sequences that contain plasmid-like replication and partition proteins as putative P-Ps.

3.  Calculate weighted gene repertoire relatedness (wGRR) between putative P-Ps and include reference dataset (e.g. typed P-Ps). Assign new P-Ps to well-defined groups, diverse communities or singletons using the similarity to the reference set P-Ps.

```{bash eval=FALSE, include=FALSE}
MM-GRC/
├── DESCRIPTION_MM-GRC.md     # describes how to use MM-GRC for P-P detection
├── Tree_based_classification.R      # R script for RF classification, predicts phage scores for plasmid-annotated elements
├── wGRR_MGEs.R                      # R script to compute wGRR between genomes
├── par_plasmid_annotation.hmm       # HMM profiles for plasmid partition proteins
├── rep_plasmid_annotation.hmm       # HMM profiles for plasmid replication proteins
├── 200122_func_cat_phage_hmms.xlsx  # Functional annotation table for phage-specific HMM profiles (phage.hmm)
├── 0523_plasmids_for_MM_GRC.xlsx    # List of plasmids used as MM-GRC input
├── model_prediction_stats_for_0523_plasmids.tsv    # Output with P-P scores from Tree_based_classification.R 
```

# Usage

## **Setup & input data**

To run MM-GRC ensure the following components are installed and accessible:

**`MM-GRC/` with necessary input data.**

**Furter relevant input data is available on [<https://zenodo.org/uploads/add_valid_link>]:**

-   `phage.hmm` – phage-specific HMM profiles.

-   `models/` – folder containing pre-trained RF models.

**External tools:**

-   [**HMMER**](<https://github.com/EddyRivasLab/hmmer>)  installed (used for hmmsearch).

-   [MMseqs2](https://github.com/soedinglab/MMseqs2) (for all-vs-all protein sequence comparison).

-   R (version ≥ 4.0) with required packages:

    -   `ranger` for random forest model and `foreach` for iterative processing
    -   `rhmmer` to parse and work with HMMER outputs and `seqinr` for biological sequence data analysis
    -   `data.table`, `readxl`, `tidyverse` for data handling and manipulation.

**Input dataset:**

-   **A protein fasta file of target sequence** for the hmmsearch. Protein sequences should be named as given by default from prodigal: genome_contig-name_protein-number (e.g. contig1_1, genome_contig1_2, genome_contig1_3, etc).

## **To run MM-GRC**

### **A) For plasmid sequences**

#### 1. Search for phage proteins

Use `hmmsearch` (from HMMER) to scan target plasmid protein sequences (e.g. `plasmid_proteins.faa`) against HMM profiles in `phage.hmm` to detect phage-specific features. 

#### 2. Detect putative P-P using the random forest models

Run the `Tree_based_classification.R` script 

***Important***: Before running, update all file paths in the script to match your local directory structure.

-   Required input:

    1.  `phage_hmm_plasmid_proteins_output.tbl.out` – HMMER output from step 1.
    2.  `200122_func_cat_phage_hmms.xlsx` – `phage.hmm` functional categories of phage HMMs.
    3.  `plasmid_proteins.faa` – protein sequences of the target plasmids.
    4.  `models/` – random forest models.

-   Output files:

    1.  `model_features_predicted_stats.txt` includes mean phage probability per plasmid (averaged across all 10 models).
    2.  `PP_list_stats_by_model.tsv` - putative P-Ps (phage probability \> 0.5).
    3.  `PP_protein_seqs.faa` - protein fasta file of predicted P-Ps (used as input for the next step).

#### 3. Run MMseqs2 to compare protein sequences

Concatenate predicted P-P protein sequences (`PP_protein_seqs.faa`) with protein sequences from a reference P-P dataset. MMseqs2 `search` will use this file (e.g., `PP_protein_seqs_extended.faa`) to perform an all-vs-all protein sequence comparison. Output file from `convertalis` (e.g., `pp_protein_seqs.m8`) is processed in the next step.

#### 4. Compute wGRR values

Run `wGRR_MGEs.R` to compute the wGRR for all P-Ps (putative and reference).

***Important***:

1.  Before running, update all file paths in the script to match your local directory structure.
2.  Set the working directory to the folder containing `wGRR_MGEs.R` script.

You can run `wGRR_MGEs.R` interactively (in RStudio) or from the command line using:

```{bash}
R < wGRR_MGEs.R --no-save
```

-   Input files:

    1.  `pp_protein_seqs.m8` — pairwise comparison output table of all P-P proteins.
    2.  `PP_protein_seqs_extended.faa` — fasta file with all P-P proteins (putative and reference).

-   Output files:

    1.  `pp_BBH_df.tsv` — best bidirectional hits (BBHs).
    2.  `pp_wGRR_df.tsv` — wGRR.

#### 5. Assigning P-P types

P-Ps are classified into P-P types based on wGRR scores.

Criteria:
1.  wGRR ≥ 0.5 to one of the known P-Ps
2.  At least 50% of its proteins match a typed P-P

Note: If a P-P matches multiple typed P-Ps, the type of the best wGRR (highest value) is assigned. 

### **B) MM-GRC for phages**

#### 1. Annotating phagen genomes with plasmid functions

Use `hmmsearch` (HMMER) is used to scan target phage protein sequences (e.g. `phages.faa`) against plasmid HMMs `rep_plasmid_annotation.hmm` for replication and `par_plasmid_annotation.hmm` for partition systems.

#### 2. Detect putative P-P

If a phage sequences contains replication and/or partition systems, it is considered as a P-P. We applied following criteria: **≥ 2 partition (par) genes**, or contain **≥ 1 replication (rep) gene and ≥ 1 partition (par) gene,** or the phage is showing a **wGRR ≥ 0.4** with known P-Ps.
We strongly recommend a manual inspection to validate these candidates.

#### Steps 3 to 5 are the same as above (computing wGRR and assign types).

# Results summary on 05/23 dataset

-   As input, we used 38,051 sequences annotated as plasmids (referred to as the “05/23” dataset). The list of NCBI accession numbers is available in `MM-GRC/0523_plasmids_for_MM_GRC.xlsx`.

-   MM-GRC was run using 25 CPUs (on the Migale computational cluster) with total runtime of ~25 h (of 23h37m were needed for the HMM search against `phage.hmm`)

-   Detection summary for 05/23 dataset (see Table S1):

    -   Total number of P-Ps new in 05/23 (and not already present in 03/21): 926

    -   Number of these P-Ps assigned to well-defined types: 498

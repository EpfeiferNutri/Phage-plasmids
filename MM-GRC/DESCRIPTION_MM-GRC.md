# Introduction

Multi-model gene repertoire clustering, or short **MM-GRC,** is a modular computational pipeline designed to detect and type phage–plasmids (P-Ps) across large phage and plasmid datasets. It can detect P-Ps from well-defined types, diverse communities and novel elements (singletons).

This document describes how to run **MM-GRC** on plasmid and phage sequences and how to interpret the output.

**MM-GRC workflow overview:**

1.  Screening of target proteins using specific HMM profiles.

2.  Initial P-P detection:

    A)  Detect P-Ps in plasmid sequences using a pre-trained random forest (RF) classifiers.
    B)  Curate phage sequences manually for the presence of plasmid-like replication and partition proteins.

3.  Calculate weighted gene repertoire relatedness (wGRR) values among P-Ps.

4.  Assign new P-Ps to well-defined groups, diverse communities or singletons based on similarity to a reference set of known P-Ps .

Note: consider elements with genome sizes between 10 kb and 300 kb.

**Folder structure:**

```{bash eval=FALSE, include=FALSE}
MM-GRC/
├── DESCRIPTION_MM-GRC.md     # describes how to use MM-GRC for P-P detection
├── Tree_based_classification.R      # R script for RF classification, predicts phage scores for plasmid-annotated elements
├── wGRR_MGEs.R                      # R script to compute wGRR between genomes
├── par_plasmid_annotation.hmm       # HMM profiles for plasmid partition proteins
├── rep_plasmid_annotation.hmm       # HMM profiles for plasmid replication proteins
├── 200122_func_cat_pvogs_pfams_tigfams.xlsx  # Functional annotation table for phage-specific HMM profiles (phage.hmm)
├── 0523_plasmids_for_MM_GRC.xlsx    # List of plasmids used as MM-GRC input
├── model_prediction_stats_for_0523_plasmids.tsv    # Output with P-P scores from Tree_based_classification.R 
```

# Usage

## **Setup & input data**

To run MM-GRC successfully, ensure the following components are installed and accessible:

**`MM-GRC/` folder with necessary input data.**

**Other input data available on [<https://zenodo.org/uploads/add_valid_link>]:**

-   `phage.hmm` – phage-specific HMM profiles.

-   `models/` – folder containing pre-trained RF models.

**External tools:**

-   HMMER installed (used for hmmsearch) [link]

-   MMseqs2 (for all-vs-all protein sequence comparison) [link]

-   R (version ≥ 4.0) with required packages:

    -   `ranger` for RF model implementation and `foreach` for iterative processing
    -   `rhmmer` to parse and work with HMMER outputs and `seqinr` for biological sequence data analysis
    -   `data.table`, `dtplyr`, `readxl`, `foreach`, `tidyverse` for data handling and manipulation.

**Input dataset:**

-   **A FASTA-formatted amino acid file of target genomes**, sutable for the hmmsearch. Protein sequences should be named as given by default from prodigal: genome/contig-name_protein-number (e.g. contig1_1, genome/contig1_2, genome/contig1_3, etc).

## **Running MM-GRC pipeline**

### **A) MM-GRC for plasmid sequences**

#### 1. Search for phage proteins in plasmid genomes

Use `hmmsearch` (from HMMER) to scan target plasmid protein sequences (e.g. `plasmid_proteins.fasta`) against HMM profiles in `phage.hmm`. to detect phage-specific features. 

#### 2. Detect putative P-P using the RF classifiers

Run the `Tree_based_classification.R` script to identify P-P candidates on phage-specific HMM hits from the previous step.

***Important***: Before running, update all file paths in the script to match your local directory structure.

-   Required input:

    1.  `phage_hmm_plasmid_proteins_output.tbl.out` – HMMER output from step 1.
    2.  `200122_func_cat_pvogs_pfams_tigfams.xlsx` – `phage.hmm` functional annotation data.
    3.  `plasmid_proteins.fasta` – protein sequences of the target plasmids.
    4.  `models/` – folder containing pre-trained random forest models.

-   Output files:

    1.  `model_features_predicted_stats.txt` includes mean phage probability per plasmid (averaged across all models)
    2.  `PP_list_stats_by_model.tsv` - filtered list of putative P-Ps (phage probability \> 0.5)
    3.  `PP_protein_seqs.fasta` - multi fasta file with protein sequences of predicted P-Ps (used as input for the next step)

#### 3. Run MMseqs2 to compare protein sequences

To assess protein-level similarity, first combine the predicted P-P protein sequences (`PP_protein_seqs.fasta`) with protein sequences from the reference P-P dataset. This combined file (e.g., `PP_protein_seqs_extended.fasta`) will be used as input for MMseqs2 `search` to perform all-vs-all protein sequence comparison and detect homologous proteins across putative and known P-Ps. Output file from `convertalis` (e.g., `pp_protein_seqs.m8`) is processed in the next step.

#### 4. Compute wGRR values

Run the `wGRR_MGEs.R` script to calculate pairwise gene repertoire relatedness for all candidate and reference P-Ps.

***Important***:

1.  Before running, update all file paths in the script to match your local directory structure.
2.  Set the working directory to the folder containing `wGRR_MGEs.R` script.

You can run `wGRR_MGEs.R` interactively (in RStudio) or from the command line using:

```{bash}
R < wGRR_MGEs.R --no-save
```

-   Input files:

    1.  `pp_protein_seqs.m8` — protein similarity output from MMseqs2.
    2.  `PP_protein_seqs_extended.fasta` — merged protein dataset containing putative and reference P-P sequences.

-   Output files:

    1.  `pp_BBH_df.tsv` — table of best bidirectional hits (BBHs).
    2.  `pp_wGRR_df.tsv` — table of pairwise weighted gene repertoire relatedness (wGRR) values.

#### 5. Assigning P-P types

Newly predicted P-Ps are classified into defined P-P types based on wGRR scores.

A candidate genome should be assigned to a known P-P type if:

1.  It has wGRR ≥ 0.5 to one of the known P-Ps
2.  At least 50% of its proteins match a known P-P.

If multiple hits are found, the P-P is assigned the type of the element with the highest wGRR value.

### **B) MM-GRC for phages**

#### 1. Search for plasmid-associated proteins in phage genomes

Use `hmmsearch` (HMMER) to scan target phage protein sequences (e.g. `phages.fasta`) against plasmis-specific profile HMMs `rep_plasmid_annotation.hmm` for replication and `par_plasmid_annotation.hmm` for partition systems.

#### 2. Detect putative P-P

Phage genomes can be considered candidate P-Ps if they meet one of the following criteria: contain **≥ 2 partition (par) genes**, or contain **≥ 1 replication (rep) gene and ≥ 1 partition (par) gene,** or the phage is showing a **wGRR ≥ 0.4** with known P-Ps

***Important***: Manual curation is required to validate these candidates.

#### Steps 3 to 5 are the same as for plasmids.

Proceed with MMseqs2 comparison, wGRR computation, and type assignment as described in "**A) MM-GRC for plasmids"** above.

# Results summary on 05/23 dataset

-   As input, we used 38,051 plasmid sequences retrieved from the non-redundant NCBI RefSeq database in May 2023 (referred to as the “05/23” dataset). All sequences were annotated as plasmids in RefSeq. The list of NCBI accession numbers is available in `MM-GRC/0523_plasmids_for_MM_GRC.xlsx`.

-   MM-GRC was run on the Migale computational cluster using 25 CPUs with total runtime \~25 hours (including 23h37m for HMM search against `phage.hmm`)

-   Detection summary for 05/23 dataset (see Table S1):

    -   Total number of new P-Ps detected: 926

    -   Number of new P-Ps assigned to well-defined P-P types: 500

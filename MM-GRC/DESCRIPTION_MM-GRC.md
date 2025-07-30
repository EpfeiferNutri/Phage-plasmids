# Introduction

Multi-model gene repertoire clustering, or short **MM-GRC,** is a modular computational pipeline developed to detect and type phage–plasmids (P-Ps) in large-scale phage and plasmid datasets.

The pipeline uses a two-part strategy:

-   plasmids are screened for phage-related signatures using phage HMM profiles and a trained random forest classifiers (RF models)

-   phages are scanned for plasmid-associated functions such as replication and partition systems.

Candidate P-Ps are further clustered based on their genome similarity to known P-Ps into well-defined groups, diverse communities or remain singletons using weighted gene repertoire relatedness (wGRR) values. For more detailed information see [REF].

This document describes the procedure for using **MM-GRC** on plasmid and phage sequences and explains how to interpret the results.

## Overview of the MM-GRC workflow:

The standard MM-GRC workflow involves the following steps:

1.  Run `hmmsearch` to annotate proteins of target sequences [REF]:

    -   use plasmid-related HMMs for phage datasets
    -   use phage-related HMMs for plasmid datasets

2.  Identify putative P-Ps using one/both of the two approaches (depending on the dataset):

    A)  Classify plasmid sequences using a pre-trained RF models (`Tree_based_classification.R`).

    B)  Curate phage sequences manually for the presence of plasmid-like replication and partition proteins.

*Note*: consider elements with genome sizes between 10 kb and 300 kb.

3.  Cluster protein sequences of all putative and reference P-Ps using `MMseqs2` [REF].

4.  Compute genome similarity using weighted gene repertoire relatedness (wGRR) scoring (`wGRR_MGEs.R`).

5.  Assign P-P types based on similarity to a reference set of known P-Ps.

# Usage

## **Setup & Installation**

To run MM-GRC successfully, ensure the following components are installed and accessible:

**External tools**:

-   HMMER (for `hmmsearch`) [link]

-   MMseqs2 (for protein sequence clustering) [link]

-   R (version ≥ 4.0) with required packages (see R script header for dependencies) [link]

**`MM-GRC/` folder contains necessary data:**

-   Required HMM files:

    -   `phage.hmm` with functional annotation table `200122_func_cat_pvogs_pfams_tigfams_v4.xlsx`

    -   `rep_plasmid_annotation.hmm`

    -   `par_plasmid_annotation.hmm`

-   R scripts `Tree_based_classification.R` and pre-trained RF models from folder `models/` for plasmid classification.

-   R scripts for genome similarity analysis (`wGRR_MGEs.R`).

-   Reference P-P list with curated type assignments.

## Input data **requirements**

-   **A FASTA-formatted amino acid file of target genomes**, sutable for the hmmsearch. Protein sequences should be named as given by default from prodigal: genome/contig-name_protein-number (e.g. contig1_1, genome/contig1_2, genome/contig1_3, etc).

## **Running MM-GRC pipeline**

### **A) MM-GRC for plasmids**

This section describes how to execute the MM-GRC pipeline on plasmid sequences.

#### 1. Search for phage proteins in plasmid genomes

Use `hmmsearch` (HMMER) to scan target plasmid protein sequences (e.g. `plasmid_proteins.fasta`) against phage-specific profile HMMs `phage.hmm`. For more detailed documentation on HMMER see [REF].

***Important:*** `phage.hmm` can be downloaded from [<https://zenodo.org/uploads/add_valid_link>]

Example command:

```{bash}
hmmsearch -o tmp.all.out --domtblout phage_hmm_plasmid_proteins_output.tbl.out phage.hmm plasmid_proteins.fasta
```

#### 2. Detect putative P-P using the RF classifiers

Run the `Tree_based_classification.R` script to identify P-P candidates on phage-specific HMM hits from the previous step.

***Important***: Before running, update all file paths in the script to match your local directory structure.

-   Required input:

    1.  `phage_hmm_plasmid_proteins_output.tbl.out` – HMMER output

    2.  `200122_func_cat_pvogs_pfams_tigfams_v4.xlsx` – `phage.hmm` functional annotations

    3.  `plasmid_proteins.fasta` – protein sequences of the target plasmids

    4.  `models/` – folder containing pre-trained random forest models

-   Output files:

    1.  `model_features_predicted_stats.txt` includes mean phage probability per plasmid (averaged across all models)
    2.  `PP_list_stats_by_model.tsv` - filtered list of putative P-Ps (phage probability \> 0.5)
    3.  `PP_protein_seqs.fasta` - multi fasta file with protein sequences of predicted P-Ps (used as input for the next step)

#### 3. Protein clustering using MMseqs2

To assess protein-level similarity, first combine the predicted P-P protein sequences (`PP_protein_seqs.fasta`) with protein sequences from the reference P-P dataset. This combined file (e.g., `PP_protein_seqs_extended.fasta`) will be used as input for MMseqs2 to perform all-vs-all clustering and detect homologous proteins across putative and known P-Ps.

Example commands:

```{bash}
mmseqs createdb  PP_protein_seqs_extended.fasta qDB/qDB
mmseqs createdb  PP_protein_seqs_extended.fasta tDB/tDB
mmseqs createindex tDB/tDB tmp -s 7
mmseqs search qDB/qDB tDB/tDB rDB/pp_protein_seqs.out tmp -s 7 -a --max-seqs 1000
mmseqs convertalis qDB/qDB tDB/tDB rDB/pp_protein_seqs.out pp_protein_seqs.m8
```

#### 4. Compute gene repertoire similarity (wGRR)

Run the `wGRR_MGEs.R` script to calculate pairwise gene repertoire relatedness for all candidate and reference P-Ps.

***Important***:

1.  Before running, update all file paths in the script to match your local directory structure.

2.  Set the working directory to the folder containing `wGRR_MGEs.R` script.

You can run `wGRR_MGEs.R` interactively or from the command line using:

```{bash}
R < wGRR_MGEs.R --no-save
```

-   Input files:

    1.  `pp_protein_seqs.m8` — protein similarity output from MMseqs2 clustering

    2.  `PP_protein_seqs_extended.fasta` — merged protein dataset containing putative and reference P-P sequences

-   Output files:

    1.  `pp_BBH_df.tsv` — table of best bidirectional hits (BBHs)

    2.  `pp_wGRR_df.tsv` — table of pairwise weighted gene repertoire relatedness (wGRR) values

#### 5. Assigning P-P types

Newly predicted P-Ps are classified into defined P-P types based on wGRR scores.

A candidate genome should be assigned to a known P-P type if:

1.  It has wGRR ≥ 0.5 to one of the known P-Ps

2.  At least 50% of its proteins match a known P-P.

If multiple hits are found, the P-P is assigned the type of the element with the highest wGRR value.

### **B) MM-GRC for phages**

This section describes how to execute the MM-GRC pipeline on phage sequences.

#### 1. Search for plasmid-associated proteins in phage genomes

Use `hmmsearch` (HMMER) to scan target phage protein sequences (e.g. `phages.fasta`) against plasmis-specific profile HMMs `rep_plasmid_annotation.hmm` for replication and `par_plasmid_annotation.hmm` for partition systems. For more detailed documentation on HMMER see [REF].

Example command:

```{bash}
hmmsearch -o tmp.all.out --domtblout rep_plasmid_hmm_phages_prt_out.tbl.out rep_plasmid_annotation.hmm phages.fasta
hmmsearch -o tmp.all.out --domtblout par_plasmid_hmm_phages_prt_out.tbl.out par_plasmid_annotation.hmm phages.fasta
```

#### 2. Detect putative P-P

Phage genomes are considered candidate P-Ps if they meet one of the following criteria:

-   contain **≥ 2 partition (par) genes**, or

-   contain **≥ 1 replication (rep) gene and ≥ 1 partition (par) gene**

-   phages showing a **wGRR ≥ 0.4** with known P-Ps

***Important***: Manual curation is required to validate these candidates.

#### Steps 3 to 5 are the same as for plasmids

Proceed with MMseqs2 clustering, wGRR computation, and type assignment as described in the "**A) MM-GRC for plasmids"** above.

## Example run: MM-GRC on 05/23 plasmid dataset

### Dataset overview

The analysis was performed on a collection of 38,051 plasmid sequences retrieved from the NCBI RefSeq database in May 2023. This dataset is referred to as the "05/23" plasmid set. All sequences were annotated as plasmids in RefSeq.

The list of NCBI accession numbers used in this study is available in the following file:

```{r}
plasmids_0523_for_MM_GRC = read_excel("MM-GRC/0523_plasmids_for_MM_GRC.xlsx")
plasmids_0523_for_MM_GRC
```

### Pipeline runnig time (approximate)

**Hardware:** 25 CPUs

Runtime **hmmsearch**: 23h 37min

Runtime **`Tree_based_classification.R`**: 26 min

Runtime **MMseqs2**: 8 min

Runtime **`wGRR_MGEs.R`**: 32 min

**Total runtime**: \~25h

### Detection summary for 05/23 plasmids

-   Total detected new P-Ps: xxxx

-   Detected from well-defined P-Ps: xxxx

## Interpretation and recommendations

-   Clustering using protein-sequence similarity such as wGRR tends to accumulate too long or too short sequences.

-   MM-GRC can detect P-Ps from well-defined types, diverse communities and novel elements (singletons)

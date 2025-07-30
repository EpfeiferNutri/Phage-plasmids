# Efficient detection and typing of phage-plasmids [REF]

This GitHub repository provides the documentation, scripts, and data necessary for the **detection of phage-plasmids (P-Ps)**. The repository is organized into several folders, each corresponding to a specific analytical approach. It also includes the data and script required to reproduce the figures and analyses presented in the associated study [REF].

For detailed information on each section, please go to the corresponding folder, which includes a **`DESCRIPTION_X.Rmd`** file describing its contents and usage.

## Repository structure:

### *tyPPing*

tyPPing is a user-friendly pipeline designed for the sensitive and efficient detection of prevalent phage-plasmid (P-P) types. It uses distinct protein profiles that were specifically generated and trained on well-characterized P-P. It classify P-Ps into defined types, along with confidence-level assessments, fast and accuratly.

Folder structure:

-   `DESCRIPTION_tyPPing.Rmd` describes how to use tyPPing for P-P detection.

-   `tyPPing.R` – main script for analyzing complete genome sequences.

-   `tyPPing_for_draft_genomes.R` – adapted script for analyzing draft genomes represented as sets of contigs.

-   `tyPPing_input_data/` contains required input data to run tyPPing:

    Supporting information tables used by tyPPing to detect and type P-Ps.

-   `small_example_test_data/` includes example input and output files for 20 test genomes processed by `tyPPing.R`.

### *MM-GRC*

MM-GRC ([m]{.underline}ulti-[m]{.underline}odel [g]{.underline}ene [r]{.underline}epertoire [c]{.underline}lustering) is an integrated approach developed to detect and classify P-Ps [REF]. This method combines functional annotation using phage- and plasmid-specific HMM profiles with machine learning (Random Forest, RF) models for P-P detection. It further employs gene repertoire relatedness to cluster and type P-Ps. It can detect all P-P types including diverse communities and novel elements.

Folder structure:

-   `DESCRIPTION_MM-GRC.Rmd` describes how to use MM-GRC for P-P detection.

-   Plasmid-specific HMMs (`par_plasmid_annotation.hmm` and `rep_plasmid_annotation.hmm`)

-   Phage-specific HMM functional annotation reference table (`200122_func_cat_pvogs_pfams_tigfams_v4.xlsx`)

-   `Tree_based_classification.R` - custom R script to compute genome annotation profiles and apply the classification models to predict phage scores ranging from 0 (plasmid-like) to 1 (phage-like) for target genomes.

-   `wGRR_MGEs.R` - custom R script used to calculate *weighted gene repertoire relatedness (wGRR)* between genomes, used for clustering and typing P-Ps.

-   P-P prediction data:

    `0523_plasmids_for_MM_GRC.xlsx` - input list of plasmids analyzed.

    `model_prediction_stats_for_0523_plasmids.tsv` - output of P-P classification from the RF models.

### *geNomad_for_P-Ps*

geNomad is a computational tool designed to classify nucleotide sequences as phages, integrated prophages, or plasmids [REF]. In this study, geNomad is applied to a set of plasmid sequences to identify potential P-Ps. This section describes the usage and interpretation of geNomad’s results in the context of P-P detection.

Folder structure:

-   `DESCRIPTION_geNomad.Rmd` describes how to use geNomad for P-P detection.

-   `0523_plasmids_for_geNomad.xlsx` - list of plasmid sequences used as input for geNomad analysis.

-   P-P detection results - geNomad output summary tables:

    `geNomad_plasmid_summary_on_0523_plasmids.tsv` - elements predicted to be plasmids

    `geNomad_virus_summary_on_0523_plasmids.tsv` - elements predicted to be phages or prophages

### *vConTACT2_for_P-Ps*

vConTACT v2 [REF] is a viral classification tool that organizes viral genomes into clusters based on gene-sharing networks. In this project, we apply vConTACT v2 to classify P-Ps by integrating query plasmid genomes (predicted to be phages by geNomad) with a curated reference dataset of known P-Ps [REF]. This part describes the usage and interpretation of vConTACT’s results in the context of P-P clustering.

Folder structure:

-   `DESCRIPTION_vConTACT2.Rmd` describes how to use vConTACT v2 for P-P detection.

-   `0523_plasmids_ref_PPs_for_vConTACT.csv` - list of genomes used as input for vConTACT v2, including target elements and reference P-Ps.

-   Clustering results (raw output from vConTACT v2):

    `genome_by_genome_overview.csv` - summary of pairwise genome associations.

    `c1.ntw` - network file representing gene-sharing relationships.

-   `reference_0321_PP_list_for_vConTACT2.xlsx` - reference P-P dataset from [REF] with assigned P-P types

-   `vConTACT2_PP_list_in_0523_plasmids.xlsx` - final output: P-Ps identified within the plasmid dataset, assigned with predicted P-P types based on clustering with the reference set.

### *Study_associated_data*

This folder contains the scripts, data, and supplementary materials needed to reproduce the figures and analyses presented in the associated publication [REF].

Folder structure:

-   `R_scripts_to_generate_figures/`

-   `Supplementary_data/` - contains all required tables

-   Documents with supplemental figures

## Big files are placed in a zenodo repository:

**[Link to the repository]**

-   `tyPPing_signature_profiles.hmm` – HMM profiles specific to P-P types, used for protein-to-profile searches by tyPPing.

-   `phage.hmm` - phage-specific HMM profiles used by MM-GRC

-   `models/` - folder contains random forest (RF) classification models trained to detect P-Ps in plasmid datasets used by MM-GRC

-   `g2g_plot_tables/` - necessary data for script `all_g_to_g_plots_filtered.R` from *Study_associated_data* folder

-   `tyPPing_criteria_tables/` - necessary data for script `figures_methods.R` from *Study_associated_data* folder

## References

-   [tyPPing]
-   [geNomad documentation](https://portal.nersc.gov/genomad/index.html)
-   [vConTACTv2 documentation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/)
-   [mm-grc]

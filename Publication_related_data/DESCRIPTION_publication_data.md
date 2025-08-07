# Reproducibility Materials for [REF]

This folder contains the scripts, data, and supplementary materials required to reproduce the figures and analyses from the publication **[REF]**.

### Folder Structure

**`R_scripts_to_generate_figures/`** folder contains all R scripts used to generate the figures presented in the publication.

**Supplemental Figures**:

-   **`supplemental document S1_tyPPing_criteria.pptx`** – figures showing the curation of cp32-related P–Ps, and testing of cutoffs for MinProteins and Composition branches and size ranges.

-   **`supplemental document S2_genome_to_genome.pptx`** – genome-to-genome plots, produced with `gggenomes` (<https://github.com/thackl/gggenomes>). Includes examples of high/medium/low confidence predictions, elements with atypical sizes, false positives/negatives.

**`Supplementary_data/`** folder contains the following data tables:

1.  **`T1_Characteristics_of_the_10_well-defined_P-P_types.xlsx`** table summarizes key features of ten well-characterized P-P types.

    -   **P-P type** - identifier for the phage-plasmid type.
    -   **Key example** - a representative, experimentally confirmed P-P element used as a reference.
    -   **Virion morphology (of example)** - virion morphology for the key example.
    -   **Taxa of host (most frequent)** - the most common host in which this P-P type is found.
    -   **Number of protein families (highly conserved)** - the total umber of highly conserved protein families defined for this P–P type. Protein sequences of these families were used to generate signature HMM profiles.
    -   **Cutoff: MinProteins** - minimum number of conserved signature proteins required for a genome to qualify as a member of this P–P type (used in the MinProteins branch of tyPPing).
    -   **Cutoff: Composition + tolerance** - expected range of composition sizes, with a gap tolerance (±X) indicating how many profiles can be missing while still considered a valid match (used in the Composition branch of tyPPing).
    -   **Size range (kb)** - typical genome size range (in kilobases) observed for P–Ps of this type, used by tyPPing to determine confidence levels.

2.  **`T2_Criteria_for_confidence_levels_used_by_tyPPing.xlsx`** table defines how the tyPPing tool assigns confidence levels (High, Medium, Low, or Not a P-P) to predicted phage–plasmid candidates based on a set of criteria.

    -   **Confidence** – final confidence level assigned to the genome.
    -   **MinProteins** – indicates whether the genome passes the MinProteins threshold, i.e., contains the minimum number of conserved signature genes.
    -   **Composition** – indicates whether the genome matches one of the expected composition of the Composition branch (with allowed tolerance for missing genes).
    -   **Match to size range** – whether the genome size falls within the expected size range for the predicted P–P type. A ±10% tolerance is allowed for borderline cases.
    -   **Remark** - Explanation of what the confidence level implies.

3.  **`T3_tyPPing_predictions_for_draft_genomes.xlsx`**

    **[To be added]**

4.  **`S1_All_MGE_information_table.xlsx`** table contains key information about all MGEs (including plasmids, phages, and putative P–Ps) used in the study. The MGEs were originally collected from non-redundant NCBI RefSeq databases (Phage_DB or Plasmid_DB) in March 2021 ("03/21" dataset) and in May 2023 ("05/23" dataset). Classification was determined using MM-GRC, geNomad, vConTACT2, and tyPPing.

    -   **NCBI accession** - accession number for the MGE sequence.
    -   **in 03/21** - indicates if the MGE was present in 03/21 dataset (Yes/No).
    -   **in 05/23** - indicates if the MGE was present in 05/23 dataset (Yes/No).
    -   **Source** - the original database source (Phage_DB or Plasmid_DB).
    -   **Host genus** - genus of the bacterial host. Host assignment was based on GenBank metadata (`ORGANISM` field) for plasmids, and the virus–host DB ([virus-host database](https://www.genome.jp/virushostdb/) for phages.
    -   **Name** - full name of the MGE.
    -   **Replicon size (bp)** - total size in base pairs.
    -   **MGE type MM-GRC** - MGE classification (e.g., Phage, Plasmid, or P–P with assigned type) predicted by the MM-GRC pipeline. For elements present in both datasets, the newer (05/23) classification is used.
    -   **MGE type geNomad** - MGE type assigned by geNomad (Phage/Prophage or Plasmid, if confidently assingned). Applied only to plasmid-dataset elements from both 03/21 and 05/23 datasets.
    -   **P‑P type vConTACT2** - P–P type based on vConTACT2 clustering. Applied to plasmid-dataset elements identified as phages by geNomad. Applied only to new elements from 05/23 since elements from 03/21 were used as references.
    -   **P‑P type tyPPing** - P–P type assigned by tyPPing. For elements present in both datasets, the newer (05/23) classification is used.
    -   **Confidence level (tyPPing)** - confidence level of the tyPPing prediction: High, Medium, or Low.
    -   **Signature to all genes ratio** - proportion of conserved signature genes (detected by tyPPing) relative to the total number of genes in the element. Applies even to elements not predicted as P–Ps but showing partial similarity.
    -   **P‑P type for the ratio** - the P–P type associated with the reported signature-to-gene ratio. If multiple types are detected, the one with the highest ratio is shown.

5.  **`S2_Profile_information_table.xlsx`** table provides detailed information on the HMM profiles defined for each P–P type, based on conserved protein families.

    -   **P–P type** - the P–P type to which the HMM profile belongs (one of the 10 well-defined types).
    -   **HMM profile ID** - identifier for the HMM profile (includes type and repsesentative protein ID).
    -   **Protein representative ID** - the representative protein sequence (as defined by MMseqs2) used to build the profile (RefSeq locus tag format).
    -   **P–P score** - score from 0 to 1 indicating how specifically this profile differentiates P–Ps from phages and plasmids (1 = highly specific).
    -   **P–P type score** - score from 0 to 1 reflecting specificity to the assigned P–P type.
    -   **Sequence score threshold** - HMMER score threshold used to determine a significant match for the MinProteins branch of tyPPing (as it is typically done for protein families in PFAMs).
    -   **Neff value** - number of effective sequences used to build the HMM; reflects diversity in the alignment (higher = more diverse).
    -   **geNomad-related columns**: geNomad hit, geNomad annotation, pident (geNomad), evalue (geNomad), bitscore (geNomad) - represents functional annotation and alignment metrics based on the geNomad database.
    -   **PHROGs-related columns**: PHROGs hit, PHROGs Annotation, Protein category (PHROGs), evalue (PHROGs), cov_profile (PHROGs)- functional annotations based on PHROGs database, including found PHROGs hit, spesific and more general function.
    -   **Plasmid specific proteins**: Rep/Par hit, evalue (Rep/Par), cov_profile (Rep/Par) - detected replication/partition proteins associated with plasmids.
    -   **Proteins for HMM (Locus tags)** - list of all protein locus tags used to build the HMM for this profile, separated by "\|\|"

6.  **tyPPing output tables** for analized datasets (for detailed column description see documentation `tyPPing/DESCRIPTION_tyPPing.md`):

    For elements from 03/21 dataset:

    -   **`All_hmm_hits_table_0321_by_tyPPing.tsv`** - contains all HMM hits detected by tyPPing, including those that do not meet the confidence criteria.
    -   **`Final_prediction_table_0321_by_tyPPing.tsv`** - lists final predictions made by tyPPing, along with assigned confidence levels (High, Medium or Low).

    For elements from 05/23 dataset:

    -   **`All_hmm_hits_table_0523_by_tyPPing.tsv`**
    -   **`Final_prediction_table_0523_by_tyPPing.tsv`**

### Other supplementary data are stored in Zenodo repository [LINK] due to storage limits:

-   **`g2g_plot_tables/`** folder contains standard input tables used by the `all_g_to_g_plots_filtered.R` script to generate genome-to-genome visualizations with **gggenomes** ([gggenomes GitHub](https://github.com/thackl/gggenomes)):

    -   **`genomes_used_gggenomes.tsv`** – contig/chromosome sequences.
    -   **`protein_info_used_gggenomes.tsv`** - protein positions and annotations.
    -   **`BBH_links_used_gggenomes.tsv`** - bidirectional best-hit (BBH) links between proteins from different genomes wich are used to connect two locations between two different sequences.
    -   **`wGRR_info_used_gggenomes.tsv`** - information about pairwise wGRR values.

-   **`tyPPing_criteria_tables/`** folder provides input data for the `figures_methods.R` script:

    -   **`all_pers_hmm_phages_plasmids_0321.tbl.out**` and** `protein_to_genome_0321.tsv`\*\* – HMM search output and corresponding protein-to-genome mapping file; used for defining MinProteins thresholds.
    -   **`cp32_0312_wGRR_df.tsv`** – wGRR values of all non-chromosomal Borrelia and Borreliella sequences (of 03/21), that are homologous to reference cp32.
    -   **`Table_all_composition_tolerance_to_gaps_cutoff.tsv`**, **`Table_P1_1_composition_pattern_size_cutoff.tsv`**, **`Table_P1_1_composition_positive_hit_cutoff.tsv`** - these tables show how different parameter settings of Composition branch affect the detection of P–P types. The table is used to assess how reliably P–P type-specific compositions can be used under varying parameter settings, particularly different tolerance to gaps.

-   **`draft_genomes_analysis/`** folder includes 12 selected draft genomes from carbapenem-resistant *Enterobacteriales* (CRE) species and predictions generated by `tyPPing_for_draft_genomes.R` for these genomes. They and organized into subfolders for short-read, long-read assemblies and complete P–P genomes.

    Each genome is provided as a nucleotide `.fasta` file. Corresponding tyPPing predictions are available in **`Final_prediction_table_XXX.tsv`**.

    A summary of these analyses is provided in: **`T3_tyPPing_predictions_for_draft_genomes.xlsx`** (described above)

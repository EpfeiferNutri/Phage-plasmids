# **tyPPing** and other methods for detecting phage-plasmids

**For more information, see: "Efficient detection and typing of phage-plasmids"** [ADD LINK]

Here you can find the documentation, scripts, and data necessary for the **identification of phage-plasmids (P-Ps)**:

<img width="2059" height="682" alt="image" src="https://github.com/user-attachments/assets/8bf089c7-8f3c-4b96-b533-60152637b8b3" />

See DESCRIPTION.md in the separate folders for details. Data and scripts required to reproduce analysis and figures of our study are in 'Publication_related_data'.

### *tyPPing*

tyPPing is a user-friendly, fast and accurate method to detect P-Ps. Currently it finds P-Ps of the type AB_1, P1_1, P1_2, N15, SSU5, pMT1, pCAV, pSLy3, pKpn, and cp32.
It uses protein profiles to search sequences for patterns (frequency and compositional sets) of conserved P-P proteins. If a match also fits the typical size range, it is predicted as a P-P with a distinct confidence.

### *MM-GRC*

MM-GRC (**m**ulti-**m**odel **g**ene **r**epertoire **c**lustering) is our first method to classify P-Ps (see PMID: 33590101). It is an integrated approach combining functional annotation (with phage- and plasmid-specific HMM profiles), machine learning (random forest) models, and which was complemented with an exhaustive literature review. It relies on the gene repertoire relatedness to type P-Ps, and detects various types including diverse communities and unrelated putative P-Ps (singletons).

### *geNomad_vConTACT2*

Here we describe how to use geNomad and vConTACT v2 for detecting and typing P-Ps.

geNomad classifies nucleotide sequences as phages, integrated prophages, or plasmids. In our study, we used it to analyze sequences of a plasmid database, cases classed as phages were considered as potential P-Ps.

vConTACT v2 is clusters viral genomes with a reference dataset using their shared gene content. Here, we used it to group the putative P-Ps identified by geNomad with 1416 P-Ps that we typed in previous work.

### *Publication_related_data*

Scripts, data, and supplementary materials to reproduce the figures and analyses presented in the publication.

### *Zenodo repository*

Further files are available at the Zenodo repository: **10.5281/zenodo.16616313**

-   `tyPPing_signature_profiles.hmm` – 763 HMM profiles (concatenated) specific to the 10 P-P types. Used for protein-to-profile comparison (needed for tyPPing).

-   `phage.hmm` – phage-specific HMM profiles (required for MM-GRC).

-   `models/` – random forest models trained to detect P-Ps in plasmid datasets (used by MM-GRC).

-   `g2g_plot_tables/` and `tyPPing_criteria_tables/` contain the tables required to produce figures of our study using the scripts `all_g_to_g_plots_filtered.R` and `figures_methods.R` (in `Publication_related_data/`).

-   `draft_genomes_analysis/` – 12 complete P-P genomes that we detected in 9 draft genomes of carbapenem-resistant *Enterobacteriales* species. We used a modified version of tyPPing, `tyPPing_for_draft_genomes.R`, on drafts assembled by short and long reads, and on hybrid assemblies.

## References

-   [MM-GRC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7969092/)
-   [HMMER](https://github.com/EddyRivasLab/hmmer)
-   [MMseqs2](https://github.com/soedinglab/MMseqs2)
-   [geNomad](https://portal.nersc.gov/genomad/index.html)
-   [vConTACTv2](https://bitbucket.org/MAVERICLab/vcontact2/src/master/)
-   [PHROGs](https://phrogs.lmge.uca.fr/)
-   [Migale](https://migale.inrae.fr/)

## Citing tyPPing
If you like tyPPing and you use it for your work, please cite: 
'Efficient detection and typing of phage-plasmids.'
[ADD LINK]



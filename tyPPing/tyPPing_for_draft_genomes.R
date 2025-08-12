##########################################################
# tyPPing: Bioinformatic Tool for P-P detection and typing 
##########################################################


suppressMessages({
  suppressWarnings({
    
    
    cat("\n-------------------  Welcome to tyPPing  ----------------------\n")
    cat("-----------------------------START-----------------------------\n\n")
    
    start.time = Sys.time()
    
    ################################################################################
    
    # Load required libraries
    
    load_libraries = function() {
      packages = c("data.table", "argparse", "tidyverse")
      for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          message(paste0("Package '", pkg, "' not found. Attempting to install..."))
          tryCatch({install.packages(pkg)}, error = function(e) {
            message(paste0("Failed to install '", pkg, "': ", e$message))})
        }
        # Try loading the package again
        if (!require(pkg, character.only = TRUE)) {
          stop(paste0("Package '", pkg, "' is not installed or failed to load. Please install it manually."))}}}
    
    load_libraries()
    
    ################################################################################
    
    # Set the required file paths
    
    if (interactive()) {
      
      ################################################################################
      # START USER INPUT SECTION
      
      # Please provide the full file paths to your own input files.
      # Edit ONLY the paths below according to your data location.
      
      ###  Input files required for your dataset analysis:
      
      # Path to the file mapping protein IDs to genome IDs (TSV format)
      protein_to_genome_path = "/path/to/your/data/protein_to_genome.tsv"
      
      # Path to the file containing genome sizes (TSV format)
      contig_size_path = "/path/to/your/data/contig_sizes.tsv"
      
      # Path to the HMMER search output file (tbl.out format)
      HMM_search_output_path = "/path/to/your/data/HMM_search_output.tbl.out"
      
      # Path to the directory where output files will be saved (should exist prior to running the script)
      output_dir_path = "/path/to/your/output/directory/"
      
      ### tyPPing parameters data path (provided with the scripts in folder `tyPPing input data`)
      
      Compositions_unique_path = '/path/to/tyPPing_parameters/Compositions_information_table.tsv'
      Profile_info_path = '/path/to/tyPPing_parameters/Profile_information_table.tsv'

      # END USER INPUT SECTION
      ################################################################################
      
    } else {
      
     # Parse command-line arguments
      parser = ArgumentParser(description = 'tyPPing: Predict P-P types from HMM search output')
      
      ###  User input files required for your dataset analysis:
      parser$add_argument("-m", "--map", required = TRUE, 
                          help = "Path to the protein-to-genome mapping file (TSV)")
      parser$add_argument("-s", "--sizes", required = TRUE, 
                          help = "Path to the genome size file (TSV)")
      parser$add_argument("-i", "--hmm_domtbl", required = TRUE, 
                          help = "Path to the HMMER output file (*.tbl.out)")
      parser$add_argument("-o", "--outdir", required = TRUE, 
                          help = "Path to the output directory")
      
      ### tyPPing parameters data path 
      parser$add_argument("--compositions", default = "tyPPing_input_data/Compositions_information_table.tsv", 
                          help = "Path to the composition information table (provided with the scripts in folder `tyPPing input data`)")
      parser$add_argument("--profiles", default = "tyPPing_input_data/Profile_information_table.tsv", 
                          help = "Path to the profile information table (provided with the scripts in folder `tyPPing input data`)")
      
      args = parser$parse_args()
      
      # Echo input paths to terminal 
      cat("Protein to genome path: ", args$map, "\n")
      cat("Genome size path:       ", args$sizes, "\n")
      cat("HMM search output path: ", args$hmm_domtbl, "\n")
      cat("Output directory path:  ", args$outdir, "\n")
      
      # Assign to working variables
      protein_to_genome_path = args$map
      contig_size_path = args$sizes
      HMM_search_output_path = args$hmm_domtbl
      output_dir_path = args$outdir
      
      Compositions_unique_path = args$compositions
      Profile_info_path = args$profiles
      
    }
    
    ################################################################################
    
    # Read and process input data
    
    cat("Reading and processing the input data...\n")
    
    # Read protein-to-genome assignments
    protein_to_genome = fread(protein_to_genome_path)
    n_protein_per_contig = protein_to_genome %>% group_by(contig_id) %>% summarise(n_protein=n())
    
    # Read genome sizes and merge with protein counts
    contig_size = fread(contig_size_path) %>% left_join(n_protein_per_contig, by="contig_id") %>% replace(is.na(.),0)
    contig_size_filtered = contig_size %>% filter(size <= 300000 & n_protein >=3)
    protein_to_genome_filtered = protein_to_genome %>% inner_join(contig_size_filtered)
    contig_to_genome_filtered = protein_to_genome_filtered %>% group_by(contig_id, genome_id) %>% summarise()
    
    # Read P-P type-specific Composition data
    Compositions_unique = fread(Compositions_unique_path) %>% 
      mutate(PP_category = case_when(`P-P type` == "AB_1" ~ "AB", `P-P type` == "P1_1" ~ "P11", `P-P type` == "P1_2" ~ "P12", T ~ `P-P type`), 
             Composition = `Composition ID`, 
             hmm_profile = `HMM profile in composition`,
             composition_size = `Composition size`)
    
    # Read signature profile information end extract necessary scores
    Profile_info = fread(Profile_info_path) 
    MinProteins_score = Profile_info %>% 
      select(hmm_profile = `HMM profile ID`, profile_scores_1 = `P-P score`, profile_scores_2 = `P-P type score`, `P-P type`) %>% 
      mutate (PP_category = case_when(`P-P type` == "AB_1" ~ "AB",  `P-P type` == "P1_1" ~ "P11", `P-P type` == "P1_2" ~ "P12", T ~ `P-P type`))
    MinProteins_cutoff = Profile_info %>% 
      select(hmm_profile = `HMM profile ID`, sequence_score_thr = `Sequence score threshold`)

    ################################################################################
    # Protein-to-profile HMM output pre-processing
    
    cat("Reading and processing the protein-to-profile HMM output file...\n")
    
    # the colnames 
    HMM_domtbl_out_colnames = c("domain_name","domain_accession","domain_len","query_name","query_accession","qlen","sequence_evalue","sequence_score",
                                "sequence_bias","domain_N","domain_of","domain_cevalue","domain_ievalue","domain_score","domain_bias"
                                ,"hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc","description")  
    
    # when reading in the table:  
    HMM_search_output = read_table(HMM_search_output_path, comment = "#", col_names = HMM_domtbl_out_colnames) %>%
      mutate(ali_length = abs(ali_to - ali_from+1), cov_profile = round(ali_length/qlen,3)) %>% 
      select(protein_id = domain_name, hmm_profile = query_name, domain_ievalue, cov_profile, sequence_score) %>%
      inner_join(protein_to_genome_filtered %>% select(protein_id, contig_id, genome_id), by="protein_id") %>%
      # count the hits per genome_id
      group_by(genome_id, hmm_profile) %>% 
      arrange(domain_ievalue) %>% 
      mutate(n_hits_profile = 1:n()) %>%
      #just take the best hit (lowest domain_ievalue)
      filter(n_hits_profile == 1) %>%
      mutate(PP_category = str_extract(hmm_profile, pattern=".*(?=_pers_)")) 
    
    cat("Protein-to-profile HMM output file processed successfully\n")
    
    ################################################################################
    # Cutoffs information table
    
    ALL_cutoffs = data.frame(
      check.names = FALSE,
      `P-P type` = c("AB_1", "cp32", "N15", "P1_1", "P1_2", "pCAV", "pKpn", "pMT1", "pSLy3", "SSU5_pHCM2"),
      composition_tol_thr = c(7, 6, 9, 11, 17, 3, 2, 6, 3, 3),
      MinProteins_min_N = c(55, 18, 7, 49, 46, 30, 38, 50, 47, 71),
      size_min = c(101329, 27653, 39838, 72057, 79071, 102639, 98807, 85066, 85190, 94107),
      size_max = c(122444, 32756, 65948, 125381, 103576, 118316, 125007, 115126, 133911, 125751),
      `10% mean size (bp)` = c(11189, 3020, 5289, 9872, 9132, 11048, 11191, 10010, 10955, 10993),
      PP_category = c("AB", "cp32", "N15", "P11", "P12", "pCAV", "pKpn", "pMT1", "pSLy3", "SSU5_pHCM2"))
    
    ################################################################################
    # MinProteins branch function: 
    # Description: 
    #   Searches for highly conserved sequences in each genome.
    #   Keeps only genomes passing a minimal threshold of required protein hits.
    
    MinProteins_score_based_prediction_f = function(category, HMM_filtered) {
      
      # keep only hits >= MinProteins cutoff values 
      HMM_filtered_MinProteins = left_join(HMM_filtered, MinProteins_cutoff %>% select(hmm_profile, sequence_score_thr), by="hmm_profile") %>% 
        filter(sequence_score >= sequence_score_thr)
      
      if (nrow(HMM_filtered_MinProteins) > 0) {
        
        # minimum required number of protein hits 
        category_MinProteins_cutoff = ALL_cutoffs$MinProteins_min_N[ALL_cutoffs$PP_category == category]
        
        # profile presence/absence matrix per genome
        check_profiles = HMM_filtered_MinProteins %>% ungroup() %>%
          select(hmm_profile, contig_id)  %>% mutate(values = 1)  %>% 
          pivot_wider(names_from = contig_id, values_from = values)  %>%  
          column_to_rownames("hmm_profile") %>% replace(is.na(.),0)
        
        profile_presence_genome_id = MinProteins_score  %>% filter(PP_category == category) %>%  select(hmm_profile, profile_scores_1, profile_scores_2) %>% 
          full_join(check_profiles %>% mutate(hmm_profile=rownames(check_profiles)), by='hmm_profile') %>% replace(is.na(.),0)
        
        profile_scores_1_genome_id_ = profile_presence_genome_id %>% select(-hmm_profile, -profile_scores_1, -profile_scores_2) * profile_presence_genome_id$profile_scores_1
        profile_scores_1_genome_id = bind_cols(MinProteins_score %>% filter(PP_category == category) %>%  select(hmm_profile, profile_scores_1, profile_scores_2), profile_scores_1_genome_id_)
        
        profile_scores_2_genome_id_ = profile_presence_genome_id %>% select(-hmm_profile, -profile_scores_1, -profile_scores_2) * profile_presence_genome_id$profile_scores_2
        profile_scores_2_genome_id = bind_cols(MinProteins_score %>% filter(PP_category == category) %>%  select(hmm_profile, profile_scores_1, profile_scores_2), profile_scores_2_genome_id_)
        
        
        TEST_info_by_genome_id_category_before_filtering = data.frame(n_hits_MinProteins = colSums(profile_presence_genome_id %>% select(-hmm_profile, -profile_scores_1, -profile_scores_2))) %>% 
          mutate(data.frame(score_1_MinProteins = colSums(profile_scores_1_genome_id %>% select(-hmm_profile, -profile_scores_1, -profile_scores_2)))) %>% 
          mutate(data.frame(score_2_MinProteins = colSums(profile_scores_2_genome_id %>% select(-hmm_profile, -profile_scores_1, -profile_scores_2)))) %>% 
          mutate(contig_id = rownames(.)) %>% left_join(contig_size, by="contig_id") %>% mutate(PP_category=category, cutoff_MinProteins = category_MinProteins_cutoff) %>%
          mutate(size = as.numeric(size), n_protein =as.numeric(n_protein), n_profiles_MinProteins= as.numeric(n_hits_MinProteins))  %>% 
          mutate(score_0_MinProteins = n_profiles_MinProteins / n_protein) %>% 
          mutate(score_0_MinProteins = round(score_0_MinProteins, digits = 3), score_1_MinProteins = round(score_1_MinProteins, digits = 1), score_2_MinProteins = round(score_2_MinProteins, digits = 1)) %>%
          select(contig_id, size, n_protein, PP_category, cutoff_MinProteins,  n_hits_MinProteins, score_0_MinProteins, score_1_MinProteins, score_2_MinProteins)
        
        TEST_info_by_genome_id_category = TEST_info_by_genome_id_category_before_filtering %>% filter(score_0_MinProteins >= 0.1)
        
        
      } else { 
        
        TEST_info_by_genome_id_category = data.frame(contig_id = "NA", size = 0, n_protein = 0, PP_category = "NA", cutoff_MinProteins = 0, n_hits_MinProteins = 0, score_0_MinProteins = 0, score_1_MinProteins = 0, score_2_MinProteins = 0)[-1, ]
      } 
      
      TEST_info_by_genome_id_category_summary = TEST_info_by_genome_id_category %>% 
        left_join(contig_to_genome_filtered) %>% group_by(genome_id, PP_category, cutoff_MinProteins) %>%
        arrange(desc(n_hits_MinProteins)) %>% 
        summarise(
          n_hits_MinProteins_list = paste(n_hits_MinProteins, collapse = ";"), 
          contig_id_list = paste(contig_id, collapse = ";"),  
          n_protein = sum(n_protein),
          size = sum(size),
          n_hits_MinProteins = sum(n_hits_MinProteins), 
          score_0_MinProteins = round(n_hits_MinProteins / n_protein, digits = 3),
          score_1_MinProteins = round(sum(score_1_MinProteins), digits = 1), 
          score_2_MinProteins = round(sum(score_2_MinProteins), digits = 1)) %>% 
        select(genome_id, size_MinProteins = size, n_protein_MinProteins = n_protein, PP_category, cutoff_MinProteins, n_hits_MinProteins, n_hits_MinProteins_list, 
               contig_id_list_MinProteins = contig_id_list, score_0_MinProteins, score_1_MinProteins, score_2_MinProteins)
      
      
      return(TEST_info_by_genome_id_category_summary)
      
    }
    
    ################################################################################ 
    
    # Composition branch function: 
    # Description: 
    #   Selects genomes matching expected P-P genome organization based on known compositions.
    #   Requires at least 50% coverage on HMM profiles & tolerance for gaps.
    
    Composition_search_based_prediction_f = function(category, HMM_filtered) {
      
      tolerance_threshold = ALL_cutoffs$composition_tol_thr[ALL_cutoffs$PP_category == category]
      
      contigs_filtered_by_hits = HMM_filtered %>% filter(cov_profile >= 0.5) %>% 
        group_by(contig_id) %>% mutate(number_hits = n()) %>% ungroup() %>%
        left_join(contig_size %>% select(-genome_id), by="contig_id") %>% mutate(score_0_c = number_hits / n_protein) %>% 
        filter(score_0_c >= 0.1) %>% select(-size, -score_0_c, -number_hits, -n_protein)
      
      composition_comparison = Compositions_unique %>% inner_join(contigs_filtered_by_hits, by = "hmm_profile") %>% 
        group_by(Composition, genome_id) %>% mutate(number_hits = n()) %>% ungroup() %>%
        filter(number_hits == composition_size | (number_hits >= composition_size - tolerance_threshold & number_hits >= 0.75*composition_size)) %>%
        group_by(contig_id) %>% arrange(desc(number_hits)) %>% slice(1) %>% ungroup() %>% 
        left_join(contig_size %>% select(-genome_id), by="contig_id") %>% arrange(desc(size)) 
      
      composition_comparison_ = composition_comparison %>% group_by(genome_id) %>% 
        mutate(contig_id_list = paste(contig_id, collapse = ";"), total_size = sum(size), total_n_protein = sum(n_protein)) %>% slice(1) %>% ungroup()%>% 
        mutate(tol_thr = tolerance_threshold)
      
      TEST_info_by_genome_id_category = composition_comparison_ %>% 
        select(genome_id, size_c = total_size, n_protein_c = total_n_protein, PP_category = PP_category.x, Composition, 
               composition_size, tol_thr, n_hits_composition = number_hits, contig_id_list_p = contig_id_list)
      
      return(TEST_info_by_genome_id_category)
      
    }
    
    ################################################################################
    # MAIN ANALYSIS LOOP: Run prediction branches by P-P type 
    
    cat("P-P prediction for each P-P type has started...\n")
    
    # Extract P-P categories present in HMM search output
    PP_category_list = HMM_search_output$PP_category %>% unique()
    
    # Initialize empty dataframes to store cumulative results
    SUMMARY_MinProteins_score_based_prediction_ALL = data.frame(NULL)
    SUMMARY_Composition_search_based_prediction_ALL = data.frame(NULL)
    
    # Loop over each P-P category and run both branches (MinProteins + Composition)
    for (category in PP_category_list) { 
      
      mapped_category = case_when(
        category == "AB"  ~ "AB_1",
        category == "P11" ~ "P1_1",
        category == "P12" ~ "P1_2",
        category == "SSU5" ~ "SSU5_pHCM2",
        TRUE              ~ category)
      
      cat("Searching for ", mapped_category, "\n")
      
      # Filter HMM hits specific for the current category
      HMM_filtered = HMM_search_output %>% filter(PP_category == category)  
      if (category == "SSU5") {
        HMM_filtered = HMM_filtered %>% mutate(PP_category = "SSU5_pHCM2")
        category = "SSU5_pHCM2"
      }
      
      # Apply MinProteins analytical branch
      MinProteins_result = MinProteins_score_based_prediction_f(category, HMM_filtered)
      # Apply Composition analytical branch
      Composition_result = Composition_search_based_prediction_f(category, HMM_filtered)
      
      # Combine results generated by P-P type into the global summary tables
      SUMMARY_MinProteins_score_based_prediction_ALL = SUMMARY_MinProteins_score_based_prediction_ALL %>% bind_rows(MinProteins_result)
      SUMMARY_Composition_search_based_prediction_ALL = SUMMARY_Composition_search_based_prediction_ALL %>% bind_rows(Composition_result)
      
    }
    
    cat("P-P prediction by P-P type finished successfully\n")
    
    
    ################################################################################
    # FINAL PREDICTION TABLE GENERATION
    
    cat("Summarizing the tyPPing predictions...\n")
    
    # Merge MinProteins and Composition results into one combined summary table
    FINAL_SUMMARY_TABLE_ = full_join(SUMMARY_MinProteins_score_based_prediction_ALL, 
                                     SUMMARY_Composition_search_based_prediction_ALL, by = c("genome_id", "PP_category")) %>% 
      mutate( PP_category = case_when(
        PP_category == "AB"  ~ "AB_1",
        PP_category == "P11" ~ "P1_1",
        PP_category == "P12" ~ "P1_2",
        PP_category == "SSU5" ~ "SSU5_pHCM2",
        TRUE              ~ PP_category))
    
    FINAL_SUMMARY_TABLE = FINAL_SUMMARY_TABLE_ %>%
      rename(`Genome ID` = genome_id, `P-P type` = PP_category, `Genome size MinProteins (bp)` = size_MinProteins, `Number of proteins MinProteins` = n_protein_MinProteins,  `Tolerance to gaps` = tol_thr,
             `MinProteins cutoff` = cutoff_MinProteins, `MinProteins hits` = n_hits_MinProteins,  `MinProteins hits list`  = n_hits_MinProteins_list, `MinProteins contigs list` = contig_id_list_MinProteins,
             `Composition size` = composition_size, `Composition hits` = n_hits_composition, `Genome size Composition (bp)` = size_c,`Composition contigs list` = contig_id_list_p, `Number of proteins Composition` = n_protein_c, 
             `Signature to all genes ratio` =  score_0_MinProteins, `P-P score sum` = score_1_MinProteins, `P-P type score sum` = score_2_MinProteins)
    
    ### Filter for predicted candidates: positive by any branch
    FINAL_PREDICTION_TABLE_ = FINAL_SUMMARY_TABLE_ %>% filter(n_hits_MinProteins >= cutoff_MinProteins | !is.na(Composition)) 

    # If any candidates found
    if (nrow(FINAL_PREDICTION_TABLE_) > 0) { 
      
      # Assign which branch predicted each result
      FINAL_PREDICTION_TABLE = FINAL_PREDICTION_TABLE_ %>%
        
        group_by(genome_id) %>% mutate(n_predictions = 1:n()) %>% 
        mutate(Predicted_by = case_when(n_hits_MinProteins >= cutoff_MinProteins & !is.na(Composition) ~ "MinProteins + Composition",
                                        n_hits_MinProteins >= cutoff_MinProteins ~ "MinProteins",
                                        !is.na(Composition) ~ "Composition"))
      
      categories_special = c("cp32", "SSU5_pHCM2", "pKpn", "pSLy3")
      categories_other = c("AB_1", "N15", "P1_1", "P1_2", "pCAV", "pMT1")
      prediction_methods = c("MinProteins + Composition", "MinProteins")
      
      # Generate final predictions
      FINAL_PREDICTION_ = FINAL_PREDICTION_TABLE %>% 
        left_join(ALL_cutoffs %>% select(`P-P type`, size_min, size_max, `10% mean size (bp)`), by=c('PP_category'='P-P type')) %>%
        
        mutate(Confidence = case_when(size_MinProteins <= size_max & size_MinProteins >= size_min & Predicted_by == "MinProteins + Composition" ~ "High",
                                      PP_category %in% categories_special & Predicted_by == "Composition" ~ "Low", # cp32 and SSU5 super 
                                      size_MinProteins > `10% mean size (bp)` + size_max | size_MinProteins < size_min -`10% mean size (bp)` ~ "Low",
                                      T ~ "Medium")) %>%
        
        mutate(prediction_score = case_when(`Confidence` == "High" ~ 3, `Confidence` == "Medium" ~ 2, `Confidence` == "Low" ~ 1)) %>%
        mutate(type_score = case_when(`PP_category` == "pSLy3" ~ 2, `PP_category` == "pKpn" ~ 1, T ~ 3))
      
      FINAL_PREDICTION_a = FINAL_PREDICTION_ %>% filter(!PP_category %in% categories_special) 
      
      FINAL_PREDICTION_b = FINAL_PREDICTION_ %>% filter(PP_category %in% categories_special) %>%
        group_by(`genome_id`) %>% 
        arrange(desc(prediction_score), desc(type_score)) %>% slice_head(n = 1) %>% 
        ungroup()

      FINAL_PREDICTION = bind_rows(FINAL_PREDICTION_a, FINAL_PREDICTION_b) %>%
        
        mutate(Confidence = case_when(PP_category %in% categories_special & Predicted_by == "Composition" & # cp32 and SSU5 super
                                        size_MinProteins < `10% mean size (bp)` + size_max & size_MinProteins > size_min - `10% mean size (bp)` ~ "Medium",
                                      T ~ Confidence)) %>%
        
        select(genome_id, PP_category, Confidence, Predicted_by, size_MinProteins, cutoff_MinProteins, n_hits_MinProteins, n_hits_MinProteins_list, contig_id_list_MinProteins, 
               Composition, size_c, composition_size, n_hits_composition, contig_id_list_p) %>%
        rename(`Genome ID` = genome_id, `P-P type` = PP_category, `Confidence level` = Confidence, `Predicted by` = Predicted_by, `Genome size MinProteins (bp)` = size_MinProteins, 
               `MinProteins cutoff` = cutoff_MinProteins, `MinProteins hits` = n_hits_MinProteins, `MinProteins hits list`  = n_hits_MinProteins_list, `MinProteins contigs list` = contig_id_list_MinProteins,
               `Genome size Composition (bp)` = size_c, `Composition size` = composition_size, `Composition hits` = n_hits_composition, `Composition contigs list` = contig_id_list_p)
      
    } else {
      
      FINAL_PREDICTION = data.frame(`Genome ID` = "NA", `P-P type` = "NA", `Confidence level` = "NA", `Predicted by` = "NA", `Genome size MinProteins (bp)` = 0, 
                                    `MinProteins cutoff` = 0, `MinProteins hits` = 0, `MinProteins hits list`  = "NA", `MinProteins contigs list` = "NA",
                                    `Genome size Composition (bp)` = 0, `Composition size` = 0, `Composition hits` = 0, `Composition contigs list` = "NA")[-1, ]
    }
  
    
    ################################################################################
    # OUTPUT FILES
    
    cat("Writing the results...\n")
    
    # Output files: user should define output_dir_path in the input section
    
    # Define output file names
    summary_file_path = file.path(output_dir_path, "All_hmm_hits_table.tsv")
    prediction_file_path = file.path(output_dir_path, "Final_prediction_table.tsv")
    
    # Write output tables
    write_tsv(FINAL_SUMMARY_TABLE, summary_file_path)
    write_tsv(FINAL_PREDICTION, prediction_file_path)
    
    cat("Results successfully written to:\n")
    cat(paste0("- ", summary_file_path, "\n"))
    cat(paste0("- ", prediction_file_path, "\n"))
    
    ################################################################################
    # tyPPing running time
    
    running.time = difftime(Sys.time(), start.time, units = "secs")
    minutes = floor(as.numeric(running.time) / 60)
    seconds = round(as.numeric(running.time)) %% 60
    
    cat("Running time: ", minutes, "minutes", seconds, "seconds\n")
    
    cat("\n------------------------------END------------------------------\n")
    cat("---------------------------------------------------------------\n\n")
    
    
  })
})

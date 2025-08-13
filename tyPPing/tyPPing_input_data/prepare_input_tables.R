library(data.table)
library(seqinr)
library(foreach)
library(tidyverse)


### Preparing protein_to_genome.tsv input file (a multi-fasta protein file, used for the hmmsearch (e.g., proteins.faa))

# protein sequences should be named as given by default from prodigal: 
# genome/contig-name_protein-number
# e.g.  genome/contig1_1, genome/contig1_2, genome/contig1_3, etc.

input_file = "path/to/input/proteins.faa"
output_file = "path/to/genome_size.tsv"

prot_info_lst = read.fasta(seqtype = "AA", as.string = T, file=input_file)

# Get protein IDs
protein_IDs = data.frame(protein_id = unlist(lapply(prot_info_lst,getName))) %>% 
  mutate(protein_id = str_split(protein_id, pattern = "\t", simplify = T)[,1]) %>% 
  mutate(protein_id = as.character(protein_id))

# Get genome IDs
genome_IDs = protein_IDs %>% distinct() %>% 
  mutate(genome_id =  str_remove(protein_id, pattern="_[0-9]{1,6}$"))

# Write the combined results to a TSV file
write_tsv(genome_IDs, output_file)


### Preparing genome_size.tsv input file (for several input .fasta files in one folder)

input_dir = "path/to/input/nucleotide_fasta_files/"
output_file = "path/to/genome_size.tsv"

# Get the list of .fasta files in the directory
fasta_files = list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)

extract_genome_size = function(file) {
  
  # Read sequences
  sequences = read.fasta(file, as.string = FALSE)
  
  # Sum the lengths of all sequences in the file (should be 1 for complete genome)
  total_length = sum(sapply(sequences, length))
  
  # Return file name (without extension) and total length
  data.frame(
    genome_id = tools::file_path_sans_ext(basename(file)),
    size = total_length,
    stringsAsFactors = FALSE)
}

# Process each fasta file and combine results
genome_sizes = foreach(file = fasta_files, .combine = rbind) %do% {
  extract_genome_size(file)}

# Write the combined results to a TSV file
write_tsv(genome_sizes, output_file)

#####

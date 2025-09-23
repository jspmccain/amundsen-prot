### looking at frag unique peptides

## loading in the list of peptides with taxonomic assigned, including the frag genus labels.
tax_quant_data_frag <- read.csv('data/all_avail_csvs_unique_peptides_only_with_tax_frag_genus.csv')

## first need to format the peptide data
tax_quant_data_frag$peptide_formatted <- gsub("\\s*\\([^\\)]+\\)","", 
                                              tax_quant_data_frag$peptide)
tax_quant_data_frag$peptide_form <- gsub('.', 
                                         replacement = '', 
                                         x = tax_quant_data_frag$peptide_formatted, 
                                         fixed = TRUE)

## then import the in silico digests
frag_kerg_only <- read.csv('data/frag_peps_only_in_frag_kerg.csv')
frag_cyl_only <- read.csv('data/frag_peps_only_in_frag_cyl.csv')
frag_species_all <- read.csv('data/frag_peps_conserved_across_species.csv')

## adding in the specificity column to these dataframes
frag_kerg_only$specificity <- rep('kerg', nrow(frag_kerg_only))
frag_cyl_only$specificity <- rep('cyl', nrow(frag_cyl_only))
frag_species_all$specificity <- rep('genus_only', nrow(frag_species_all))

## aggregate all together
frag_pep_df <- rbind(frag_kerg_only, frag_cyl_only, frag_species_all)

## then examine how many of the peptides are expected based on the frag genomes alone 
# tax_quant_data_frag_w_specificity <- tax_quant_data_frag %>% 
#   inner_join(frag_pep_df %>% 
#               dplyr::rename(peptide_form = peptide), 
#             by = 'peptide_form')

## Evaluate how many unique, Frag. peptides there are (7077 unique Frag-specific peptides in total from the metaP data)
peptide_number_eval <- tax_quant_data_frag %>% 
  dplyr::filter(assigned_frag != 'non-frag-genus')
peptide_number_eval$peptide_form %>% unique() %>% length()

## Evaluate how many of these peptides are also present in the cultured Frag proteomic composition
peptide_to_genome_mapping <- tax_quant_data_frag %>% 
  dplyr::filter(assigned_frag != 'non-frag-genus') %>% 
  left_join(frag_pep_df %>% 
               dplyr::rename(peptide_form = peptide), 
             by = 'peptide_form')

## If the specificity column is an NA value, it was not in the set of Frag genomes used. Designate that in the specificty column:
peptide_to_genome_mapping[is.na(peptide_to_genome_mapping$specificity),]$specificity <- 'not-in-genomes'

## Evaluating how many were found in the genomes, how many map to Frag genus only, and how many map to Kerg. vs. Cyl.
peptide_to_genome_mapping %>% 
  group_by(specificity) %>% 
  summarize(number_specific_peps = n())
## 1503 cylindrus peptides that are specific
## 1807 kerguelensis peptides that are specific
## 746 peptides that are common across the genus.

## Make a smaller dataframe that has just the unique mappings
peptide_to_genome_mapping_unique_maps <- peptide_to_genome_mapping %>% 
  dplyr::select(peptide_form, specificity) %>% 
  unique()

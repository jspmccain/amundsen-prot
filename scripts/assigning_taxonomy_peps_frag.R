### ASSIGNING TAXONOMIC STRINGS TO PEPTIDES 

library(dplyr)
library(magrittr)
library(readxl)

source("scripts/post_processing_functions_amundsen_prot.R")

### reading in the taxonomic info for TARA stuff
tara_tax <- read_excel("../antiox-review/data/tara_ocean_smags/Table_S03_statistics_nr_SMAGs_METdb.xlsx", 
                       sheet = 1, skip = 2)

tara_tax$tax_string <- paste(tara_tax$Best_taxonomy_KINGDON, tara_tax$Best_taxonomy_PHYLUM, tara_tax$Best_taxonomy_CLASS, 
                             tara_tax$Best_taxonomy_ORDER, tara_tax$Best_taxonomy_FAMILY, tara_tax$Best_taxonomy_GENRE,
                             sep = ";")

### reading in the data from TFG
annot_mcl <- read_excel("../ross-sea-meta-omics/data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx",
                        sheet = 1)

all_unique_quant_peps <- read.csv("data/all_avail_csvs_unique_peptides_only.csv")

tax_string_test_frag_included <- get_tax_assoc_comprehensive_frag(quant_out = all_unique_quant_peps,
                                                                  tfg_annot = annot_mcl, 
                                                                  tara_out = tara_tax)

## Remove all spaces from the Frag. string that describes the assigned Frag taxonomy
tax_string_test_frag_included$assigned_frag <- gsub(",", "_", tax_string_test_frag_included$assigned_frag)

## appending those values to the unique peptide dataframe itself
# all_unique_quant_peps$assigned_tax_frag <- tax_string_test_frag_included$frag_genus
# all_unique_quant_peps$taxonomic_group <- tax_string_test_frag_included$taxonomic_group

write.csv(tax_string_test_frag_included, "data/all_avail_csvs_unique_peptides_only_with_tax_frag_genus.csv", row.names = FALSE, quote = FALSE)





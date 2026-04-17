### looking at taxonomic break down

library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggrepel)
library(stringr)

## peptide diatom antioxidant data
diatom_assigned_clusters <- read.csv("data/diatom_assigned_clusters_antioxidant_output.csv")
diatom_assigned_clusters_frag_only_mnsod <- diatom_assigned_clusters %>% filter(grepl("Fragilariopsis", tax_assignments),
                                                                          cluster_id == 5)
diatom_assigned_clusters_frag_only_cuznsod <- diatom_assigned_clusters %>% filter(grepl("Fragilariopsis", tax_assignments),
                                                                          cluster_id == 6)

## peptide taxonomic and functional mappings 
tax_quant_data <- read.csv("data/all_avail_csvs_unique_peptides_only_with_tax.csv")
tax_quant_data_frag <- read.csv('data/all_avail_csvs_unique_peptides_only_with_tax_frag_genus.csv')
func_quant_data <- read.csv("data/all_unique_quant_peps_w_annots_no_extended_annots.csv")

# loading peptide quant data
all_avail_csvs_normed_averaged <- read.csv("data/all_avail_csvs_normed_averaged.csv")
all_avail_csvs_normed_not_averaged <- read.csv("data/all_avail_csvs_normed_not_averaged.csv")

# making the sample ID a character
all_avail_csvs_normed_averaged$sample_id <- all_avail_csvs_normed_averaged$sample_id %>% as.character()
all_avail_csvs_normed_not_averaged$sample_id <- all_avail_csvs_normed_not_averaged$sample_id %>% as.character()

## load in the protein amounts per filter
protein_amounts <- read.csv("data/Amundsen2018_2018_BCA_all_pages.csv")
protein_amounts$sample_id <- as.character(protein_amounts$sample_id)

# reading in the metal data
sample_with_metal <- read.csv("data/protein_sample_sheets_with_metal_data_manual_selection.csv",
                              fileEncoding = "Latin1", check.names = F)
# formatting sample ID names
sample_names_alone <- str_split(string = sample_with_metal$sample_label, pattern = "_")
sample_names_alone_formatted <- unlist(lapply(sample_names_alone, `[[`, 3))
sample_names_alone_formatted[sample_names_alone_formatted == "227"] <- "227"
# sample_names_alone_formatted[sample_names_alone_formatted == "227b"] <- "227"
sample_with_metal$sample_id <- str_remove(sample_names_alone_formatted, "^0+")
sample_with_metal_subset <- sample_with_metal %>%
  dplyr::filter(sample_id %in% unique(all_avail_csvs_normed_averaged$sample_id))

### getting Korean CTD data merged:
kopri_data <- read.csv("data/kopri_data_subset.csv")
names(kopri_data) <- c("station", "depth_kopri", "po4_umol", "po4_std", "no2no3_umol", 
                       "no2no3_std", "nh4_umol", "nh4_std", "sio2_umol", "sio2_std")

merge_korean_data_with_tm_data <- function(){
  
  kopri_dataframe <- data.frame(station_kor = numeric(),
                                depth_kopri = numeric(),
                                no2no3_umol = numeric(), 
                                sio2_umol = numeric(),
                                depth_difference = numeric())
  
  for(row_i in 1:nrow(sample_with_metal)){
    # row_i <- 30
    print(row_i)
    station_nioz <- sample_with_metal[row_i, ]$station
    depth_nioz <- sample_with_metal[row_i, ]$depth_m
    station_filtered_kopri <- kopri_data %>% dplyr::filter(station == station_nioz)
    
    if(nrow(station_filtered_kopri) > 0){
      
      depth_differences <- abs(depth_nioz - station_filtered_kopri$depth_kopri)
      depth_lowest <- depth_differences == min(depth_differences)
      closest_kopri_match <- station_filtered_kopri[depth_lowest,]
      
      ## if there are two equadistant depths, take the mean value
      if(nrow(closest_kopri_match) > 1){
        closest_kopri_match <- closest_kopri_match %>% 
          group_by(station) %>% 
          summarize_all(mean)
      }
      
      # making the kopri dataframe
      kopri_dataframe_i <- data.frame(station_kor = station_nioz,
                                      depth_kopri = closest_kopri_match$depth_kopri,
                                      no2no3_umol = closest_kopri_match$no2no3_umol, 
                                      sio2_umol = closest_kopri_match$sio2_umol,
                                      depth_difference = abs(closest_kopri_match$depth_kopri - depth_nioz))
    } else if(nrow(station_filtered_kopri) == 0){
      kopri_dataframe_i <- data.frame(station_kor = NA,
                                      depth_kopri = NA,
                                      no2no3_umol = NA, 
                                      sio2_umol = NA,
                                      depth_difference = NA)
    }
    
    kopri_dataframe <- rbind(kopri_dataframe, kopri_dataframe_i)
    print(nrow(kopri_dataframe))
  }
  return(kopri_dataframe)
}

### note that the KOPRI data doesn't have station 4.
kopri_data_out <- merge_korean_data_with_tm_data()

sample_with_metal_with_korean <- cbind(sample_with_metal, kopri_data_out) %>% 
  dplyr::filter(sample_label != 'PRO_3.0_227b')

write.csv(sample_with_metal_with_korean, 
          file = 'data/intermediate_data/environmental-data-per-sample-id.csv', 
          row.names = FALSE)

sample_with_metal_with_korean[sample_with_metal_with_korean$sample_id == '227a',]$sample_id <- '227'

## merging the taxonomic, functional, and environmental data together with peptide abundances
all_avail_csvs_normed_averaged_tax_append <- all_avail_csvs_normed_not_averaged %>% 
  inner_join(tax_quant_data_frag %>% dplyr::select(-protein, -n_proteins), by = "peptide") %>% 
  inner_join(func_quant_data %>% dplyr::select(-protein, -n_proteins), by = "peptide") %>% 
  inner_join(sample_with_metal_with_korean, by = "sample_id")

### getting the proportion of the proteome that is attributed to different taxonomic groups
total_taxonomix_biomass_per_filter <- all_avail_csvs_normed_averaged_tax_append %>% 
  group_by(assigned_tax, sample_id, file_name) %>% ## file name is included here because we want an assessment of measurement error
  summarize(summed_prot_biomass = sum(norm_abundance_not_sum_to_one_dino),
            count_inf = n())

## getting normalization factors
total_taxonomix_biomass_per_filter_norm_factors <- all_avail_csvs_normed_averaged_tax_append %>% 
  group_by(sample_id, file_name) %>% ## file name is included here because we want an assessment of measurement error
  summarize(summed_prot_biomass_norm_factor = sum(norm_abundance_not_sum_to_one_dino))

total_taxonomix_biomass_per_filter_norm_factors_no_ambig <- all_avail_csvs_normed_averaged_tax_append %>% 
  filter(assigned_tax != "non-unique-grpnorm-taxgrp",
         assigned_tax != "mixture-of-taxonomic-annotations") %>% 
  group_by(sample_id, file_name) %>% ## file name is included here because we want an assessment of measurement error
  summarize(summed_prot_biomass_norm_factor_no_ambig = sum(norm_abundance_not_sum_to_one_dino)) %>% 
  ## join the two, to get the multiplication factor required
  inner_join(total_taxonomix_biomass_per_filter_norm_factors, by = c('sample_id', 'file_name')) %>% 
  dplyr::mutate(multiplication_norm_factor = summed_prot_biomass_norm_factor/summed_prot_biomass_norm_factor_no_ambig)

total_taxonomic_biomass_no_unambiguous_renorm <- all_avail_csvs_normed_averaged_tax_append %>%
  ## non assigned peptides are removed, they lead to underestimates of proportion of taxon specific biomass
  filter(assigned_tax != "non-unique-grpnorm-taxgrp",
         assigned_tax != "mixture-of-taxonomic-annotations") %>% 
  group_by(assigned_tax, sample_id, file_name) %>% 
  summarize(summed_prot_biomass_no_unknowns = sum(norm_abundance_not_sum_to_one_dino)) %>% 
  ## now attach the normalization factors, and renormalize
  inner_join(total_taxonomix_biomass_per_filter_norm_factors_no_ambig, by = c('sample_id', 'file_name')) %>% 
  dplyr::mutate(summed_prot_biomass_no_unk_renorm = summed_prot_biomass_no_unknowns*multiplication_norm_factor)

write.csv(total_taxonomic_biomass_no_unambiguous_renorm, 
          file = 'data/intermediate_data/total_taxonomic_biomass_no_unambiguous_renorm.csv', 
          quote = FALSE, row.names = FALSE)

## Checking that the renormalization worked as expected:
# tester1 <- total_taxonomic_biomass_no_unambiguous %>% 
#   group_by(sample_id, file_name) %>% 
#   summarize(test = sum(summed_prot_biomass_no_unk_renorm))
# 
# round(tester1$test, 4) == round(total_taxonomix_biomass_per_filter_norm_factors$summed_prot_biomass_norm_factor, 4)

## The challenge now is that these abundances do not sum to 1, so need to recalculate normalization factors and then renormalize.
# total_taxonomic_biomass_no_unambiguous_norm_factors <- total_taxonomic_biomass_no_unambiguous %>% 
#   ungroup() %>% 
#   group_by(sample_id, 
#            file_name) %>% 
#   summarize(norm_factors_removed_unambiguous = sum(summed_prot_biomass_no_unknowns))
# 
# ## Now rejoin these with the taxon abundance estimates from previously, and renormalize
# total_taxonomic_biomass_no_unambiguous_renorm <- total_taxonomic_biomass_no_unambiguous %>% 
#   inner_join(total_taxonomic_biomass_no_unambiguous_norm_factors, by = c('sample_id', 'file_name')) %>% 
#   dplyr::mutate(summed_prot_biomass_no_unk_renorm = summed_prot_biomass_no_unknowns/norm_factors_removed_unambiguous)

## Getting an idea of the measurement error
### First need to dcast this dataframe
replicate_dataframe_identifier <- total_taxonomic_biomass_no_unambiguous_renorm %>% 
  ungroup() %>% 
  dplyr::select(sample_id, file_name) %>% 
  unique()

write_replicates <- function(numbers_list){
  occurence_table <- table(replicate_dataframe_identifier$sample_id) %>% as.data.frame()
  replicate_vector_out <- c()
  for(i in 1:nrow(occurence_table)){
    if(occurence_table$Freq[i] == 2){
      replicate_vector_out <- c(replicate_vector_out, c("A", "B"))
    }
    if(occurence_table$Freq[i] == 1){
      replicate_vector_out <- c(replicate_vector_out, c("A"))
    }
  }
  return(replicate_vector_out)
}

replicates_to_append <- write_replicates(replicate_dataframe_identifier$sample_id)
replicate_dataframe_identifier$replicate_id <- as.factor(replicates_to_append)

total_taxonomic_biomass_no_unambiguous_across_reps <- total_taxonomic_biomass_no_unambiguous_renorm %>% 
  inner_join(replicate_dataframe_identifier, by = c('sample_id', 'file_name')) %>% 
  dplyr::select(-file_name) %>% 
  ungroup() %>%
  dcast(assigned_tax + sample_id ~ replicate_id, value.var = "summed_prot_biomass_no_unk_renorm")

## Visualizing the measurement error for Fragilariopsis abundance estimates across samples.
replicate_values_frag <- total_taxonomic_biomass_no_unambiguous_across_reps %>% 
  dplyr::filter(assigned_tax == 'Fragilariopsis') %>%
  ggplot(aes(x = A, y = B)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  xlab('Replicate A') +
  ylab('Replicate B') +
  geom_abline(slope = 1, intercept = 0) +
  # geom_label_repel(aes(label = sample_id)) +
  theme_bw();replicate_values_frag

ggsave(plot = replicate_values_frag, filename = 'figures/tester.png')

# formatting the kerg vs. cyl data ----------------------------------------

source('scripts/frag-species-genus-mapping.R')

## make a new copy of this other df. This is from 'merging peptides with quant.R' It contains the peptide intensity information
## and the sample ID information for different peptides.
all_avail_csvs_normed_averaged_with_genus_spec <- all_avail_csvs_normed 

## format the peptide sequence itself
all_avail_csvs_normed_averaged_with_genus_spec$peptide_formatted <- gsub("\\s*\\([^\\)]+\\)","", 
                                                                         all_avail_csvs_normed_averaged_with_genus_spec$peptide)
all_avail_csvs_normed_averaged_with_genus_spec$peptide_form <- gsub('.', replacement = '', 
                                                                    x = all_avail_csvs_normed_averaged_with_genus_spec$peptide_formatted, 
                                                                    fixed = TRUE)

## add in the genus_specificity stuff
all_avail_csvs_normed_averaged_with_genus_spec_2 <- all_avail_csvs_normed_averaged_with_genus_spec %>% 
  left_join(peptide_to_genome_mapping_unique_maps,
            by = 'peptide_form')

## Looking at numbers of unique peptides for cylindrus, genus only, and kerguelensis for each sample
all_avail_csvs_normed_cyl_kerg_unique_peps <- all_avail_csvs_normed_averaged_with_genus_spec_2 %>% 
  group_by(specificity, sample_id, file_name) %>% 
  summarize(sum_norm_abundance = sum(norm_abundance_not_sum_to_one_dino),
            number_of_unique_peps = n()) %>% 
  dcast(sample_id + file_name ~ specificity, value.var = 'number_of_unique_peps')

all_avail_csvs_normed_cyl_kerg_abundance <- all_avail_csvs_normed_averaged_with_genus_spec_2 %>% 
  group_by(specificity, sample_id, file_name) %>% 
  summarize(sum_norm_abundance = sum(norm_abundance_not_sum_to_one_dino),
            number_of_unique_peps = n()) %>% 
  dcast(sample_id + file_name ~ specificity, value.var = 'sum_norm_abundance') %>% 
  dplyr::mutate(fraction_kerg = kerg/(kerg + cyl))

## Join the abundance information with the number of unique peptides dataframe
all_avail_csvs_normed_cyl_kerg_abundance_w_peps <- all_avail_csvs_normed_cyl_kerg_abundance %>% 
  inner_join(all_avail_csvs_normed_cyl_kerg_unique_peps %>% 
               dplyr::rename(cyl_unique_peps = cyl,
                             genus_unique_peps = genus_only,
                             kerg_unique_peps = kerg),
             by = c('sample_id', 'file_name'))

write.csv(all_avail_csvs_normed_cyl_kerg_abundance_w_peps, 
          file = 'data/intermediate_data/all_avail_csvs_normed_cyl_kerg_abundance_w_peps.csv', 
          quote = FALSE, row.names = FALSE)

# adding in the different df's together -----------------------------------

## need to calculate the mean and SD value for each estimated biomass
total_taxonomic_biomass_no_unambiguous_renorm_summarized <- total_taxonomic_biomass_no_unambiguous_renorm %>% 
  ## inner_join with the frag biomass dataframe
  inner_join(all_avail_csvs_normed_cyl_kerg_abundance_w_peps, by = c('sample_id', 'file_name')) %>% 
  ## multiply the prot biomasses by the fraction of kerguelensis. Note this only makes the Fragilariopsis assigned_tax interpretable.
  dplyr::mutate(summed_prot_biomass_no_unk_renorm_frac = summed_prot_biomass_no_unk_renorm*fraction_kerg) %>% 
  group_by(assigned_tax, 
           sample_id) %>% 
  ## join this with the 
  summarize(mean_prot_proportion = mean(summed_prot_biomass_no_unk_renorm_frac),
            sd_prot_proportion = sd(summed_prot_biomass_no_unk_renorm_frac))

write.csv(total_taxonomic_biomass_no_unambiguous_renorm_summarized, 
          file = 'data/intermediate_data/total_taxonomix_biomass_per_filter.csv', 
          quote = FALSE, row.names = FALSE)





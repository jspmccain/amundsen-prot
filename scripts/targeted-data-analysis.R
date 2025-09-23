### targeted data analysis

library(ggplot2)
library(dplyr)
library(magrittr)
library(readxl)
library(reshape2)
library(GGally)
library(stringr)
library(ggpubr)

setwd('~/My Drive (jspmccain@gmail.com)/Projects/amundsen-prot/')

# reading in data ---------------------------------------------------------

## read in the peptide amounts
peptide_amounts <- read_excel(path = "data/peptide_amounts_left.xlsx", sheet = 4) %>% 
  dplyr::rename(sample_id = `short sample`,
                ms_sample = sample,
                ug_loaded = `ug total protein loaded`)

## read in the frag specific peptides
frag_peptides <- read_excel(path = 'data/Extra_diatom_peak_areas_20231129.xlsx') %>% 
  dplyr::rename(peak_area = `peak area`)

### read in the peptide to protein connection
peptide_to_protein <- read_excel(path = 'data/Amudsen_SRM_data_2021/Sample_SRM_export_01.xlsx', sheet = 2)

### read in the metaproteomic data that got assigned taxonomic groups
taxonomic_specific_biomass <- read.csv("data/intermediate_data/total_taxonomix_biomass_per_filter.csv")

## read in the proportion of frag species biomass per sample id
sample_id_frag_biomass <- read.csv('data/intermediate_data/all_avail_csvs_normed_cyl_kerg_abundance_w_peps.csv')

### read in environmental data associated with a single filter
env_filter_data <- read.csv("data/protein_sample_sheets_with_metal_data.csv") %>% 
  dplyr::select(time_collected, volume_filtered, sample_id, Depth,
                station, Latitude, longitude, depth_m_tm_bottle, salinity_tm_bottle,
                c.n, molar_c.n, pon_concentration_mg.ml, poc_concentration_mg.ml,
                po4_.um., si_.um., dmn_.nm., dfe_.nm., 
                par_from_protein_bottles, temp_from_protein_bottles) %>% 
  dplyr::rename(lat = Latitude,
                lon = longitude,
                depth_m_from_protein_bottles = Depth)

# format the frag peptide targeted data --------------------------------------------

## first need to parse the sample id label
str_split_list_samples <- str_split(frag_peptides$sample, "_")

## get the second last number within each sublist
sample_label_vector <- sapply(str_split_list_samples, '[[', 4)

## if the first character is a zero, remove
remove_first_char_zero <- function(list_of_sample_names){
  ### this function goes through every sample and removes a zero if the first character is zero
  
  ## make an output vector
  output_vec <- c()
  
  for(i in 1:length(list_of_sample_names)){
    sample_name_i <- list_of_sample_names[i]
    first_value <- substr(sample_name_i, start = 1, stop = 1)
    if(first_value == "0"){
      sample_name_i <- substr(sample_name_i, start = 2, stop = nchar(sample_name_i))
    }
    
    output_vec <- c(output_vec,
                    sample_name_i)
  }
  return(output_vec)
}

## append this to the frag_peptides dataframe
frag_peptides$sample_id <- remove_first_char_zero(sample_label_vector)

## format the sample loading name colum
frag_peptides$ms_sample <- substr(frag_peptides$sample, start = 10, stop = nchar(frag_peptides$sample))


# format the frag discovery data ------------------------------------------
frag_discovery_biomass <- taxonomic_specific_biomass %>% 
  dplyr::filter(assigned_tax == 'Fragilariopsis') %>% 
  dplyr::rename(summed_prot_biomass_m = mean_prot_proportion,
                summed_prot_biomass_sd = sd_prot_proportion)

## need to append the number of replicates
frag_discovery_biomass$number_reps <- ifelse(test = is.na(frag_discovery_biomass$summed_prot_biomass_sd), yes = 1, no = 2)
frag_discovery_biomass_se <- frag_discovery_biomass %>% 
  dplyr::mutate(summed_prot_biomass_se = summed_prot_biomass_sd/sqrt(number_reps))

average_standard_error <- mean(frag_discovery_biomass_se$summed_prot_biomass_se, na.rm = TRUE)

### Now add in the estimated standard error for the two samples without a replicate injection.
frag_discovery_biomass_se$summed_prot_biomass_se[frag_discovery_biomass_se$sample_id == 74] <- average_standard_error
frag_discovery_biomass_se$summed_prot_biomass_se[frag_discovery_biomass_se$sample_id == 197] <- average_standard_error

# checking overlapping sample id labels -----------------------------------

### there are three dataframes. One from the targeted (with standards), one from the targeted (without standards), one from the discovery.

## 42 samples from discovery
disc_sample_ids <- frag_discovery_biomass_se$sample_id %>% unique()

## 43 samples from targeted
targeted_st_ids <- peptide_amounts$sample_id %>% unique()

## 39 from targetd without standards
targeted_no_st_ids <- frag_peptides$sample_id %>% unique()

overlapping_sample_id <- intersect(targeted_no_st_ids, intersect(disc_sample_ids, targeted_st_ids))

## part of the issue here is that these samples sometimes have the letter B in it, which refers to a sample injection with 
## different amounts of total protein.

# look to see if the targeted sum of peptides correlates with the discovery Frag biomass estimates --------

## first need to pair the targeted Frag dataset with the targeted standard dataset. This will allow me to multiply the peak
## area by the protein loaded.

### make a dataframe that maps the sample_id to the total protein loaded
loaded_to_sample_id <- peptide_amounts %>% 
  dplyr::select(sample_id, ug_loaded) %>% 
  unique()

### joining this to frag peptides so we know how much protein is added per sample
frag_peptides_w_prot <- frag_peptides %>% 
  inner_join(loaded_to_sample_id, 
             by = "sample_id")

## adding a column where the sample ID removes the "B"
frag_peptides_w_prot$sample_id_no_b <- gsub(pattern = "B", 
                                            replacement = "", 
                                            x = frag_peptides_w_prot$sample_id)

## calculate an abundance metrics
frag_peptides_targeted_sum <- frag_peptides_w_prot %>% 
  mutate(peak_area_per_protein_loaded = peak_area/ug_loaded) %>% 
  group_by(sample_id_no_b) %>% 
  summarize(median_pappl = median(peak_area_per_protein_loaded),
            mean_pappl = mean(peak_area_per_protein_loaded),
            sd_pappl = sd(peak_area_per_protein_loaded),
            sum_pappl = sum(peak_area_per_protein_loaded),
            n_peps_detected = n()) %>% 
  dplyr::rename(sample_id = sample_id_no_b) %>% 
  mutate(sample_id = as.numeric(sample_id))

## join this targeted frag df with discovery
sum_sum_relationship_for_targeted_discovery <- frag_discovery_biomass_se %>% 
  inner_join(frag_peptides_targeted_sum, by = "sample_id") %>% 
  ggplot(aes(x = summed_prot_biomass_m, y = sum_pappl)) +
  geom_point() +
  xlab('Sum of Frag-specific peptides from Discovery MetaP') +
  ylab('Summed Peak Area per Protein Loaded') +
  geom_smooth(method = 'lm') +
  theme_bw()

mean_sd_relationship_for_targeted <- frag_discovery_biomass_se %>% 
  inner_join(frag_peptides_targeted_sum, by = "sample_id") %>% 
  ggplot(aes(x = mean_pappl, y = sd_pappl)) +
  geom_point() +
  xlab('Mean Peak Area per Protein Loaded') +
  ylab('SD Peak Area per Protein Loaded') +
  geom_smooth(method = 'lm') +
  theme_bw()

targeted_discovery_connection <- ggarrange(sum_sum_relationship_for_targeted_discovery,
                                           mean_sd_relationship_for_targeted)
ggsave(targeted_discovery_connection, filename = 'figures/targeted_discovery_connection.png')

# formatting targeted data to be normalized -------------------------------

## we need to make another column that does not have a B, because the discovery data didn't have this
peptide_amounts$sample_id_no_b <- gsub(pattern = "B", 
                                       replacement = "", 
                                       x = peptide_amounts$sample_id) %>% as.numeric()

## format the peptide amounts
peptide_amounts_form <- peptide_amounts %>% 
  dplyr::mutate(fmol_per_ug = `L/H`*`fmol heavy on column`/ug_loaded) %>% 
  inner_join(peptide_to_protein, by = 'Peptide')

## filter for antioxidant proteins
peptide_amounts_form_anti <- peptide_amounts_form %>% 
  dplyr::filter(protein %in% c('CuZnSOD_01',
                               'CuZnSOD_02',
                               'MnSOD_01',
                               'Peroxi_01',
                               'Peroxi_02',
                               'Peroxi_03',
                               'Glutathi_01'))

peptide_amounts_form_anti$protein_ag <- as.factor(peptide_amounts_form_anti$protein)  
levels(peptide_amounts_form_anti$protein_ag) <- c("cuznsod", 
                                                         "cuznsod", 
                                                         "gluta", 
                                                         "mnsod", 
                                                         "perox1", 
                                                         "perox2", 
                                                         "perox3")

# looking at protein by protein relationships -----------------------------
peptide_amounts_anti_summarized_disc <- peptide_amounts_form_anti %>% 
  dplyr::select(sample_id, protein_ag, protein, Peptide, fmol_per_ug) %>% 
  group_by(sample_id, protein_ag) %>% 
  summarize(mean_fmol_per_ug = mean(fmol_per_ug),
            sd_fmol_per_ug = sd(fmol_per_ug),
            number_reps = n()) %>% 
  dplyr::mutate(se_fmol_per_ug = sd_fmol_per_ug/sqrt(number_reps)) %>% 
  inner_join(frag_discovery_biomass_se %>% ## need to convert the sample_id
               dplyr::mutate(sample_id = as.character(sample_id)), 
             by = 'sample_id') %>% 
  dplyr::select(-assigned_tax)

## making a wide dataframe for Stan for the fmol per ug measurements
peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean <- peptide_amounts_anti_summarized_disc %>% 
  ungroup() %>% 
  dcast(sample_id ~ protein_ag, value.var = 'mean_fmol_per_ug')
peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd <- peptide_amounts_anti_summarized_disc %>% 
  ungroup() %>% 
  dcast(sample_id ~ protein_ag, value.var = 'se_fmol_per_ug')

## making a wide dataframe for Stan for the mass fraction estimates
peptide_amounts_anti_summarized_disc_wide_disc_sum_mean <- peptide_amounts_anti_summarized_disc %>% 
  ungroup() %>% 
  dcast(sample_id ~ protein_ag, value.var = 'summed_prot_biomass_m')
peptide_amounts_anti_summarized_disc_wide_disc_sum_sd <- peptide_amounts_anti_summarized_disc %>% 
  ungroup() %>% 
  dcast(sample_id ~ protein_ag, value.var = 'summed_prot_biomass_se')

sample_id_frag_biomass_more_than_50_kerg_peps <- sample_id_frag_biomass %>% 
  dplyr::filter(kerg_unique_peps > 50)

## key data outputs to be fit in Stan
# peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean
# peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd
# peptide_amounts_anti_summarized_disc_wide_disc_sum_mean
# peptide_amounts_anti_summarized_disc_wide_disc_sum_sd

## subsetting for only the high abundance of kerg
peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset <- peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean %>% 
  dplyr::filter(sample_id %in% sample_id_frag_biomass_more_than_50_kerg_peps$sample_id) %>% 
  dplyr::select(sample_id, cuznsod, mnsod, perox2)
peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset <- peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd %>% 
  dplyr::filter(sample_id %in% sample_id_frag_biomass_more_than_50_kerg_peps$sample_id) %>% 
  dplyr::select(sample_id, cuznsod, mnsod, perox2)
peptide_amounts_anti_summarized_disc_wide_disc_sum_mean_subset <- peptide_amounts_anti_summarized_disc_wide_disc_sum_mean %>% 
  dplyr::filter(sample_id %in% sample_id_frag_biomass_more_than_50_kerg_peps$sample_id) %>% 
  dplyr::select(sample_id, cuznsod, mnsod, perox2)
peptide_amounts_anti_summarized_disc_wide_disc_sum_sd_subset <- peptide_amounts_anti_summarized_disc_wide_disc_sum_sd %>% 
  dplyr::filter(sample_id %in% sample_id_frag_biomass_more_than_50_kerg_peps$sample_id) %>% 
  dplyr::select(sample_id, cuznsod, mnsod, perox2)

### append the environmental data associated with each measurement
peptides_norm_env_data <- peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean %>% 
  inner_join(env_filter_data,
             by = "sample_id")


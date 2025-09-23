library(dplyr)
library(magrittr)
library(quantreg)
library(ggrepel)
library(readxl)

# get all the global functions
source('scripts/post_processing_functions_amundsen_prot.R')

## reading in ORF annotation data:

### reading in the taxonomic info for TARA stuff
tara_tax <- read_excel("../antiox-review/data/tara_ocean_smags/Table_S03_statistics_nr_SMAGs_METdb.xlsx", 
                       sheet = 1, skip = 2)

tara_tax$tax_string <- paste(tara_tax$Best_taxonomy_KINGDON, tara_tax$Best_taxonomy_PHYLUM, tara_tax$Best_taxonomy_CLASS, 
                             tara_tax$Best_taxonomy_ORDER, tara_tax$Best_taxonomy_FAMILY, tara_tax$Best_taxonomy_GENRE,
                             sep = ";")

## reading in the data from TFG
annot_mcl <- read_excel("../ross-sea-meta-omics/data/mcmurdo-metatrans/Bertrand_McCrow_TFG/annotation_allTFG.grpnorm_mmetsp_fc_pn_reclassified.edgeR.xlsx",
                        sheet = 1)

### first I want to read in all the peptide-quant openms output

## getting the spectral count files
all_files_data <- dir('data/mzML-converted/')

# getting the spectral counts data
all_files_data_xic <- all_files_data[grepl(pattern = 'feature_xic', 
                                                   x = all_files_data)]
# reading in all the files
all_avail_csvs <- get_multiple_files(list_of_files = paste0('data/mzML-converted/', 
                                                            all_files_data_xic))

# dim(all_avail_csvs)
# names(all_avail_csvs)

## summing up intensities by file to get normalization factors
normalization_factors_by_injection <- all_avail_csvs %>% 
  group_by(file_string) %>% 
  summarize(normalization_factor_sum = sum(abundance),
            number_of_psms = n())

### read in the DB independent normalization facotrs
db_ind_norm_factors <- read.csv("data/amundsen_prot_normalization_factors.csv")
db_dino_norm_factors <- read.csv('data/dino-converted/dino_summaries.csv')

## plotting the correlation between TIC and db-dependent norm factors

# need to make consistent formatting of sample ID to merge with other data
file_names_w_mzml <- gsub(db_ind_norm_factors$file_name, 
     pattern = "../data/mzML-converted/", replacement = "", 
     fixed = TRUE)
file_names_w_o_mzml <-  gsub(file_names_w_mzml, 
                             pattern = ".mzML", replacement = "", 
                             fixed = TRUE)
db_ind_norm_factors$file_name_formatted <- file_names_w_o_mzml

# need to make consistent formatting of sample ID to merge with other data, this time with the dinosaur values
file_names_w_mzml_dino <- gsub(db_dino_norm_factors$file_name, 
                               pattern = "../data/dino-converted/", replacement = "", 
                               fixed = TRUE)
file_names_w_o_mzml_dino <-  gsub(file_names_w_mzml_dino, 
                             pattern = ".features.tsv", replacement = "", 
                             fixed = TRUE)
db_dino_norm_factors$file_name_formatted <- file_names_w_o_mzml_dino

# Now converting the openms output
db_dependent_string_w_feature <- gsub(normalization_factors_by_injection$file_string, 
                            pattern = "data/mzML-converted/", replacement = "", 
                            fixed = TRUE)
db_dependent_string_w_o_feature <- gsub(db_dependent_string_w_feature, 
                            pattern = ".featureXML_feature_xic.csv", replacement = "", 
                            fixed = TRUE)
normalization_factors_by_injection$file_name_formatted <- db_dependent_string_w_o_feature

# Merging the db indp and db dp norm factors
norm_factors_merged <- normalization_factors_by_injection %>% 
  inner_join(db_ind_norm_factors, by = "file_name_formatted") %>% 
  inner_join(db_dino_norm_factors %>% 
               dplyr::rename(file_name_tsv = file_name), by = "file_name_formatted")

## making a modelled normalization factor

## testing out quantile regression

# plot(mpg ~ wt, data = mtcars, pch = 16, main = "mpg ~ wt")
# abline(lm(mpg ~ wt, data = mtcars), col = "red", lty = 2)
# abline(rq(mpg ~ wt, data = mtcars, tau = 0.95), col = "blue", lty = 2)
# legend("topright", legend = c("lm", "rq"), col = c("red", "blue"), lty = 2)

## model the 95 quantile
rqfit_sum_pep_tic <- rq(normalization_factor_sum ~ tic, data = norm_factors_merged, tau = 0.95)

## add this to the normalization factor dataframe
norm_factors_merged$normalization_factor_sum_95 <- predict(object = rqfit_sum_pep_tic, newdata = norm_factors_merged)

# ### plotting correlation of normalization factors -----------------------

norm_factor_plot <- norm_factors_merged %>% 
  ggplot(aes(x = tic, y = normalization_factor_sum)) +
  geom_point() +
  theme_bw() +
  xlab('Total Ion Current') +
  ggtitle('A.') +
  # geom_label_repel(aes(label = file_name_formatted)) +
  ylab('Sum of Identified Peptide Intensities') +
  geom_abline(slope = 1, intercept = 0) +
  scale_size_continuous('Number of PSMs');norm_factor_plot

norm_factor_plot_dino <- norm_factors_merged %>% 
  ggplot(aes(x = tic, y = sum_of_int)) +
  geom_point() +
  theme_bw() +
  ggtitle('B.') +
  xlab('Total Ion Current') +
  geom_abline(slope = 1, intercept = 0) +
  ylab('Sum of Peptide-like Feature Intensities') +
  scale_size_continuous('Number of PSMs')

g_norm_comp <- ggarrange(norm_factor_plot, norm_factor_plot_dino)

ggsave(g_norm_comp, 
       filename = 'figures/normalization_comparison_dino_tic_sum_pep_int.svg', height = 4, width = 8)

ggsave(norm_factor_plot, filename = 'figures/norm-factor-plot.png')
ggsave(plot = norm_factor_plot_dino, filename = 'figures/norm-factor-plot-dino.png')

# appending normalization factors -----------------------------------------

all_avail_csvs_normed <- all_avail_csvs %>% 
  inner_join(norm_factors_merged, 
             by = "file_string") %>% 
  mutate(norm_abundance_pep = abundance/normalization_factor_sum,
         norm_abundance_not_sum_to_one = abundance/tic,
         norm_abundance_not_sum_to_one_dino = abundance/sum_of_int,
         norm_abundance_quantile = abundance/normalization_factor_sum_95)

## make normalization factors based on TIC-normalized abundances
tic_normalized_abundances_factors <- all_avail_csvs_normed %>% 
  group_by(file_string) %>% 
  summarize(sum_of_tic_norm = sum(norm_abundance_not_sum_to_one),
            sum_of_dino_norm = sum(norm_abundance_not_sum_to_one_dino),
            sum_of_quant_norm = sum(norm_abundance_quantile),
            sum_of_db_dep_norm = sum(norm_abundance_pep))

## now I want to merge peptides that were observed on replicate injections
all_avail_csvs_normed_averaged <- all_avail_csvs_normed %>% 
  group_by(peptide, sample_id,
           n_proteins, protein) %>% 
  summarize(mean_abundance_norm = mean(norm_abundance_not_sum_to_one_dino),
            mean_abundance_no_norm = mean(abundance),
            number_of_double_obs_peps = n())

# get tic by sample -------------------------------------------------------

extract_number <- function(file_path) {
  # Use regular expression to match the pattern "SM_" followed by two or three digits
  matches <- regmatches(file_path, regexpr("(?<=SM_)\\d{2,3}", file_path, perl = TRUE))
  
  # Return the matched value
  return(matches)
}

norm_factors_merged$sample_id_val <- extract_number(file_path = norm_factors_merged$file_string)
norm_factors_merged_mean <- norm_factors_merged %>% 
  group_by(sample_id_val) %>% 
  summarize(mean_tic = mean(tic))

write.csv(x = norm_factors_merged_mean, file = 'data/intermediate_data/mean_tic_by_sample.csv', row.names = FALSE)

# filter the peptides for only the unique ones, this is for taxonomic assignment
all_avail_csvs_unique_only <- all_avail_csvs_normed_averaged %>%
  ungroup() %>% 
  dplyr::select(peptide, 
                protein, 
                n_proteins) %>% 
  unique()

write.csv(all_avail_csvs_unique_only, "data/all_avail_csvs_unique_peptides_only.csv", row.names = FALSE, quote = FALSE)
write.csv(all_avail_csvs_normed_averaged, "data/all_avail_csvs_normed_averaged.csv", row.names = FALSE, quote = FALSE)
write.csv(all_avail_csvs_normed, "data/all_avail_csvs_normed_not_averaged.csv", row.names = FALSE, quote = FALSE)



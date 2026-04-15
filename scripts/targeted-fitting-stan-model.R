### fit the interaction model to the proteomic data
## this is the key dataframe needed from 'targeted-data-analysis.R'

# load
library(rstan)
library(tidybayes)
library(bayesplot)
library(beepr)
library(truncnorm)

# incorporating measurement error and covariates --------------------------

fit_mv_norm_proteomic_me_cov <- stan_model("scripts/cov_matrix_measurement_error_covariates_kerg.stan")

sample_number_name_df_me <- data.frame(number_val = 1:(ncol(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset) - 1),
                                       protein_name = names(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset %>%
                                                              dplyr::select(-sample_id)))
sample_index_to_sample_id <- data.frame(sample_index = 1:nrow(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset),
                                        sample_id = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$sample_id)
## for notes on this model: https://yingqijing.medium.com/multivariate-normal-distribution-and-cholesky-decomposition-in-stan-d9244b9aa623

sample_number_name_df_me$protein_name_nice <- c('CuZnSOD', 'MnSOD', 'Peroxiredoxin')

## adding in the time of day

# Convert to POSIXct
env_filter_data$time_collected_form <- as.POSIXct(env_filter_data$time_collected, 
                                                          format = "%Y-%m-%d %H:%M:%S")
env_filter_data$time_of_day <- format(env_filter_data$time_collected_form,
                                              format = "%H:%M:%S")

convert_to_hours <- function(time_strings) {
  # Convert the time strings into POSIXct format
  times <- as.POSIXct(time_strings, format = "%H:%M:%S", tz = "UTC")
  
  # Extract hours, minutes, and seconds
  hours <- as.numeric(format(times, "%H"))
  minutes <- as.numeric(format(times, "%M"))
  seconds <- as.numeric(format(times, "%S"))
  
  # Calculate the total hours since midnight
  total_hours <- hours + (minutes / 60) + (seconds / 3600)
  
  return(total_hours)
}

env_filter_data$hours_since_midnight <- convert_to_hours(env_filter_data$time_of_day)

env_filter_data_ordered <- peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset %>% 
  inner_join(env_filter_data, by = "sample_id")

## these data come from https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.00763/full#supplementary-material
liefer_data <- data.frame(pg_c_per_cell = c(0.277, 0.362, 0.582, 0.634, 15.1, 26.7, 139, 368), 
                          pg_protein_per_cell = c(0.287, 0.144, 0.573, 0.337, 13.1, 6.86, 103, 45.7)) %>% 
  dplyr::mutate(protein_per_carbon = pg_protein_per_cell/pg_c_per_cell)

## read in the Joli et al 2023 New Phyto data
joli_data <- read_excel('data/Joli-et-al2023-new-phyto-supp-values.xlsx', sheet = 2)

## Taking only the rows that have both protein and carbon data.
joli_data_cc <- joli_data[complete.cases(joli_data),] %>% 
  dplyr::mutate(protein_per_carbon = protein_ug_per_cell/(doc_ug_carbon_per_cell))

# changing prior input for metallation ------------------------------------

estimate_pars_me_cov_10_1 <- sampling(fit_mv_norm_proteomic_me_cov, 
                             data = list(N = nrow(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset), 
                                         K = 3,
                                         # P = 8,
                                         P_joli = length(joli_data_cc$protein_per_carbon),
                                         cuznsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$cuznsod,
                                         mnsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$mnsod,
                                         perox2_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$perox2,
                                         cuznsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$cuznsod,
                                         mnsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$mnsod,
                                         perox2_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$perox2,
                                         ## This looks strange -- it's on purpose, For each protein the normalization constant is the same.
                                         frag_proportion_meas = peptide_amounts_anti_summarized_disc_wide_disc_sum_mean_subset$cuznsod,
                                         frag_proportion_se = peptide_amounts_anti_summarized_disc_wide_disc_sum_sd_subset$cuznsod,
                                         temp_vals = scale(env_filter_data_ordered$temp_from_protein_bottles) %>% as.numeric(),
                                         fe_vals = scale(env_filter_data_ordered$dfe_.nm.) %>% as.numeric(),
                                         mn_vals = scale(env_filter_data_ordered$dmn_.nm.) %>% as.numeric(),
                                         light_vals = scale(env_filter_data_ordered$par_from_protein_bottles) %>% as.numeric(),
                                         # liefer_protein_per_carbon = liefer_data$protein_per_carbon,
                                         joli_protein_per_carbon = joli_data_cc$protein_per_carbon,
                                         liefer_mean_protein_per_carbon = mean(liefer_data$protein_per_carbon),
                                         liefer_sd_protein_per_carbon = sd(liefer_data$protein_per_carbon),
                                         # hours_mid = sin(env_filter_data_ordered$hours_since_midnight*pi/24),
                                         shape1 = 10,
                                         shape2 = 1), 
                             iter = 6000, 
                             chains = 4)


estimate_pars_me_cov_1_1 <- sampling(fit_mv_norm_proteomic_me_cov, 
                                     data = list(N = nrow(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset), 
                                                 K = 3,
                                                 # P = 8,
                                                 P_joli = length(joli_data_cc$protein_per_carbon),
                                                 cuznsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$cuznsod,
                                                 mnsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$mnsod,
                                                 perox2_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$perox2,
                                                 cuznsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$cuznsod,
                                                 mnsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$mnsod,
                                                 perox2_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$perox2,
                                                 ## This looks strange -- it's on purpose, For each protein the normalization constant is the same.
                                                 frag_proportion_meas = peptide_amounts_anti_summarized_disc_wide_disc_sum_mean_subset$cuznsod,
                                                 frag_proportion_se = peptide_amounts_anti_summarized_disc_wide_disc_sum_sd_subset$cuznsod,
                                                 temp_vals = scale(env_filter_data_ordered$temp_from_protein_bottles) %>% as.numeric(),
                                                 fe_vals = scale(env_filter_data_ordered$dfe_.nm.) %>% as.numeric(),
                                                 mn_vals = scale(env_filter_data_ordered$dmn_.nm.) %>% as.numeric(),
                                                 light_vals = scale(env_filter_data_ordered$par_from_protein_bottles) %>% as.numeric(),
                                                 # liefer_protein_per_carbon = liefer_data$protein_per_carbon,
                                                 joli_protein_per_carbon = joli_data_cc$protein_per_carbon,
                                                 liefer_mean_protein_per_carbon = mean(liefer_data$protein_per_carbon),
                                                 liefer_sd_protein_per_carbon = sd(liefer_data$protein_per_carbon),
                                                 # hours_mid = sin(env_filter_data_ordered$hours_since_midnight*pi/24),
                                                 shape1 = 1,
                                                 shape2 = 1), 
                                     iter = 6000, 
                                     chains = 4)

estimate_pars_me_cov_2_1 <- sampling(fit_mv_norm_proteomic_me_cov, 
                                     data = list(N = nrow(peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset), 
                                                 K = 3,
                                                 # P = 8,
                                                 P_joli = length(joli_data_cc$protein_per_carbon),
                                                 cuznsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$cuznsod,
                                                 mnsod_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$mnsod,
                                                 perox2_meas = peptide_amounts_anti_summarized_disc_wide_fmol_ug_mean_subset$perox2,
                                                 cuznsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$cuznsod,
                                                 mnsod_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$mnsod,
                                                 perox2_se = peptide_amounts_anti_summarized_disc_wide_fmol_ug_sd_subset$perox2,
                                                 ## This looks strange -- it's on purpose, For each protein the normalization constant is the same.
                                                 frag_proportion_meas = peptide_amounts_anti_summarized_disc_wide_disc_sum_mean_subset$cuznsod,
                                                 frag_proportion_se = peptide_amounts_anti_summarized_disc_wide_disc_sum_sd_subset$cuznsod,
                                                 temp_vals = scale(env_filter_data_ordered$temp_from_protein_bottles) %>% as.numeric(),
                                                 fe_vals = scale(env_filter_data_ordered$dfe_.nm.) %>% as.numeric(),
                                                 mn_vals = scale(env_filter_data_ordered$dmn_.nm.) %>% as.numeric(),
                                                 light_vals = scale(env_filter_data_ordered$par_from_protein_bottles) %>% as.numeric(),
                                                 # liefer_protein_per_carbon = liefer_data$protein_per_carbon,
                                                 joli_protein_per_carbon = joli_data_cc$protein_per_carbon,
                                                 liefer_mean_protein_per_carbon = mean(liefer_data$protein_per_carbon),
                                                 liefer_sd_protein_per_carbon = sd(liefer_data$protein_per_carbon),
                                                 # hours_mid = sin(env_filter_data_ordered$hours_since_midnight*pi/24),
                                                 shape1 = 2,
                                                 shape2 = 1), 
                                     iter = 6000, 
                                     chains = 4)

# aggregating coefficient dataframes --------------------------------------

beta_0_df <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(beta_0[n]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), by = 'n') %>% 
  dplyr::rename(beta_val = beta_0)

beta_1_df <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(beta_1_temp[n]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), by = 'n') %>% 
  dplyr::rename(beta_val = beta_1_temp)

beta_2_df <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(beta_2_fe[n]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), by = 'n') %>% 
  dplyr::rename(beta_val = beta_2_fe)

beta_3_df <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(beta_3_mn[n]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), by = 'n') %>% 
  dplyr::rename(beta_val = beta_3_mn)

beta_4_df <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(beta_4_light[n]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), by = 'n') %>% 
  dplyr::rename(beta_val = beta_4_light)

## now add a column with the coefficient identifier, and change the name
beta_0_df$coef_type <- rep('Intercept', 3)
beta_1_df$coef_type <- rep('Temperature', 3)
beta_2_df$coef_type <- rep('Iron', 3)
beta_3_df$coef_type <- rep('Manganese', 3)
beta_4_df$coef_type <- rep('Light', 3)

## aggregate all these into one coef dataframe
all_beta_df <- rbind(beta_0_df, beta_1_df, beta_2_df, beta_3_df, beta_4_df)

## getting the observed standard deviations
temp_sd <- sd(env_filter_data_ordered$temp_from_protein_bottles) %>% round(2)
dfe_sd <- sd(env_filter_data_ordered$dfe_.nm.) %>% round(2)
dmn_sd <- sd(env_filter_data_ordered$dmn_.nm.) %>% round(2)
par_sd <- sd(env_filter_data_ordered$par_from_protein_bottles) %>% round(2)

## plotting the mean value
intercept_plot <- all_beta_df %>% 
  dplyr::filter(coef_type == 'Intercept') %>% 
  ggplot(aes(x = beta_val, y = protein_name_nice)) +
  geom_point(aes(colour = protein_name_nice)) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper, colour = protein_name_nice), 
                 height = 0.2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  theme(legend.position = 'none') +
  scale_colour_manual(values = c('darkblue', 'darkorange', 'deepskyblue3')) +
  labs(y = "",
       x = expression("Intercept Coefficient (fmol Protein "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*")"));intercept_plot

ggsave(intercept_plot, filename = 'figures/intercept-antioxidant-proteins-mvnorm.pdf', height = 4.6, width = 11)

all_coefficiants_plot <- all_beta_df %>% 
  dplyr::filter(coef_type != 'Intercept') %>% 
  ggplot(aes(x = beta_val, y = protein_name_nice)) +
  geom_point(aes(colour = protein_name_nice)) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper, colour = protein_name_nice), 
                 height = 0.2) +
  facet_wrap(~coef_type, nrow = 5) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10), strip.text = element_text(size = 12)) +
  geom_vline(xintercept = 0, alpha = 0.6) +
  theme(legend.position = 'none') +
  scale_colour_manual(values = c('darkblue', 'darkorange', 'deepskyblue3')) +
  geom_label(data = data.frame(coef_type = c('Iron', 'Light', 'Manganese', 'Temperature'),
                               sd_val = c(dfe_sd, par_sd, dmn_sd, temp_sd),
                               beta_val = rep(40, 4),
                               protein_name_nice = rep('Peroxiredoxin', 4)),
            aes(label = paste('\u03C3 = ', sd_val)), alpha = 0.4) +
  labs(y = "",
       x = expression("Coefficient (fmol Protein "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*")"));all_coefficiants_plot

ggsave(all_coefficiants_plot, filename = 'figures/coefficients-antioxidant-proteins-mvnorm.svg', width = 7, height = 7)


# looking at estimated values ---------------------------------------------
estimated_norm_vals <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(y[sample_index, number_val]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me, 
             by = 'number_val')

## aggregating the estimated values mean, upper, and low estimates
mean_estimate_vals <- estimated_norm_vals %>%
  dplyr::select(-number_val) %>% 
  dcast(sample_index ~ protein_name, value.var = 'y')

lower_estimate_vals <- estimated_norm_vals %>%
  dplyr::select(-number_val) %>% 
  dcast(sample_index ~ protein_name, value.var = '.lower')
names(lower_estimate_vals) <- paste0(names(lower_estimate_vals), "_lower")

upper_estimate_vals <- estimated_norm_vals %>%
  dplyr::select(-number_val) %>% 
  dcast(sample_index ~ protein_name, value.var = '.upper')
names(upper_estimate_vals) <- paste0(names(upper_estimate_vals), "_upper")

estimated_vals_formatted <- mean_estimate_vals %>% 
  inner_join(lower_estimate_vals %>% 
               dplyr::rename(sample_index = sample_index_lower),
             by = 'sample_index') %>% 
  inner_join(upper_estimate_vals %>% 
               dplyr::rename(sample_index = sample_index_upper),
             by = "sample_index") %>% 
  inner_join(sample_index_to_sample_id, by = 'sample_index') %>% 
  inner_join(env_filter_data, by = "sample_id")


# dependence of each protein value with depth -----------------------------
mnsod_profiles <- estimated_vals_formatted %>% 
  ggplot(aes(x = mnsod, y = depth_m_from_protein_bottles)) +
  geom_point() +
  facet_wrap(~station) +
  scale_y_reverse() +
  theme_bw() +
  geom_errorbarh(aes(xmin = mnsod_lower, xmax = mnsod_upper)) +
  labs(y = "Depth (m)",
       # x = expression("fmol Fe, MnSOD "*mu*"g Fragilariopsis Protein"^-1*""),
       x = expression("fmol Fe, MnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*")"))

cuznsod_profiles <- estimated_vals_formatted %>% 
  ggplot(aes(x = cuznsod, y = depth_m_from_protein_bottles)) +
  geom_point() +
  facet_wrap(~station) +
  scale_y_reverse() +
  theme_bw() +
  geom_errorbarh(aes(xmin = cuznsod_lower, xmax = cuznsod_upper)) +
  labs(y = "Depth (m)",
       # x = expression("fmol Fe, MnSOD "*mu*"g Fragilariopsis Protein"^-1*""),
       x = expression("fmol CuZnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*")"))


peroxi_profiles <- estimated_vals_formatted %>% 
  ggplot(aes(x = perox2, y = depth_m_from_protein_bottles)) +
  geom_point() +
  facet_wrap(~station) +
  scale_y_reverse() +
  theme_bw() +
  geom_errorbarh(aes(xmin = perox2_lower, xmax = perox2_upper)) +
  labs(y = "Depth (m)",
       # x = expression("fmol Fe, MnSOD "*mu*"g Fragilariopsis Protein"^-1*""),
       x = expression("fmol Peroxiredoxin "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*")"))

profiles_separated <- ggarrange(mnsod_profiles, cuznsod_profiles, peroxi_profiles)

## some acrobatics needed to make a profile plot across stations
estimated_vals_formatted_melt_mean  <- estimated_vals_formatted %>% 
  dplyr::select(mnsod, cuznsod, perox2, depth_m_from_protein_bottles, sample_id) %>% 
  melt(c('depth_m_from_protein_bottles', 'sample_id'), variable.name = 'protein_name', value.name = 'protein_mean')

estimated_vals_formatted_melt_lower  <- estimated_vals_formatted %>% 
  dplyr::select(mnsod_lower, cuznsod_lower, perox2_lower, depth_m_from_protein_bottles, sample_id) %>% 
  melt(c('depth_m_from_protein_bottles', 'sample_id'), variable.name = 'protein_name', value.name = 'protein_lower') %>% 
  dplyr::mutate(protein_name = gsub(pattern = '_lower', replacement = '', x = protein_name)) %>% 
  dplyr::select(-depth_m_from_protein_bottles)

estimated_vals_formatted_melt_upper  <- estimated_vals_formatted %>% 
  dplyr::select(mnsod_upper, cuznsod_upper, perox2_upper, depth_m_from_protein_bottles, sample_id) %>% 
  melt(c('depth_m_from_protein_bottles', 'sample_id'), variable.name = 'protein_name', value.name = 'protein_upper') %>% 
  dplyr::mutate(protein_name = gsub(pattern = '_upper', replacement = '', x = protein_name)) %>% 
  dplyr::select(-depth_m_from_protein_bottles)

estimated_vals_formatted_melt_combined <- estimated_vals_formatted_melt_mean %>% 
  inner_join(estimated_vals_formatted_melt_lower, by = c('protein_name', 'sample_id')) %>% 
  inner_join(estimated_vals_formatted_melt_upper, by = c('protein_name', 'sample_id')) %>% 
  inner_join(env_filter_data, by = 'sample_id') %>% 
  inner_join(sample_number_name_df_me, by = 'protein_name')

## making a profile plot with all proteins on the sample plot
profiles_proteins_together <- estimated_vals_formatted_melt_combined %>% 
  dplyr::mutate(station_name = paste0("Station: ", station)) %>% 
  ggplot(aes(x = protein_mean, y = depth_m_from_protein_bottles.x)) +
  geom_point(aes(colour = protein_name_nice)) +
  facet_wrap(~station_name, nrow = 2) +
  scale_y_reverse() +
  theme_bw() +
  geom_errorbarh(aes(xmin = protein_lower, xmax = protein_upper, colour = protein_name_nice)) +
  # labs(y = "Depth (m)",
  #      x = expression("fmol Protein "*mu*"g F. kerguelensis Protein"^-1*"")) +
  labs(y = "Depth (m)",
       x = expression("fmol Protein "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +

  scale_colour_manual(values = c('darkblue', 'darkorange', 'deepskyblue3')) +
  theme(legend.title=element_blank(), legend.position = 'bottom',
        legend.margin = margin(c(0,0,-2,0)),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"), plot.margin = margin(c(-2,0,0,0)));profiles_proteins_together

ggsave(profiles_proteins_together, filename = 'figures/protein-profiles-long.pdf',
       width = 13.5*0.7, height = 4.89*0.7)

# making pairwise comparison figure ---------------------------------------

p1 <- estimated_vals_formatted %>% 
  inner_join(sample_index_to_sample_id, by = "sample_index") %>% 
  ggplot(aes_string(y = "cuznsod", x = "mnsod")) +
  geom_point(shape = 20, size = 2, alpha = 0.3) + 
  # geom_label(aes(label = sample_index)) +
  geom_errorbarh(aes_string(xmin = paste0("mnsod", "_lower"), 
                            xmax = paste0("mnsod", "_upper")),
                 alpha = 0.2, lwd = 1.05) +
  geom_errorbar(aes_string(ymin = paste0("cuznsod", "_lower"), 
                           ymax = paste0("cuznsod", "_upper")),
                alpha = 0.2, lwd = 1.05) +
  theme_bw() +
  ylab(expression("fmol CuZnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  xlab(expression("fmol MnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  geom_abline(slope = 1, 
              intercept = 0) +
  theme(text = element_text(size = 7), 
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.1,0.2,0.1,1), 'lines'))

p2 <- estimated_vals_formatted %>% 
  inner_join(sample_index_to_sample_id, by = "sample_index") %>% 
  ggplot(aes_string(y = "cuznsod", x = "perox2")) +
  geom_point(shape = 20, size = 2, alpha = 0.3) + 
  # geom_label(aes(label = sample_index)) +
  geom_errorbarh(aes_string(xmin = paste0("perox2", "_lower"), 
                            xmax = paste0("perox2", "_upper")),
                 alpha = 0.2, lwd = 1.05) +
  geom_errorbar(aes_string(ymin = paste0("cuznsod", "_lower"), 
                           ymax = paste0("cuznsod", "_upper")),
                alpha = 0.2, lwd = 1.05) +
  theme_bw() +
  ylab(expression("fmol CuZnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  xlab(expression("fmol Peroxiredoxin "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  geom_abline(slope = 1, 
              intercept = 0) +
  theme(text = element_text(size = 7), 
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.1,0.2,0.1,1), 'lines'))

p3 <- estimated_vals_formatted %>% 
  inner_join(sample_index_to_sample_id, by = "sample_index") %>% 
  ggplot(aes_string(y = "mnsod", x = "perox2")) +
  geom_point(shape = 20, size = 2, alpha = 0.3) + 
  # geom_label(aes(label = sample_index)) +
  geom_errorbarh(aes_string(xmin = paste0("perox2", "_lower"), 
                            xmax = paste0("perox2", "_upper")),
                 alpha = 0.2, lwd = 1.05) +
  geom_errorbar(aes_string(ymin = paste0("mnsod", "_lower"), 
                           ymax = paste0("mnsod", "_upper")),
                alpha = 0.2, lwd = 1.05) +
  theme_bw() +
  ylab(expression("fmol MnSOD "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  xlab(expression("fmol Peroxiredoxin "*mu*"g "*italic("F. kerguelensis")*" Protein"^-1*"")) +
  geom_abline(slope = 1, 
              intercept = 0) +
  theme(text = element_text(size = 7), 
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.1,0.2,0.1,1), 'lines'))

blank_plot <- ggplot() + theme_void()

all_by_all_est_anti <- ggarrange(p1, p2, blank_plot, p3, align = 'hv')

ggsave(all_by_all_est_anti, filename = 'figures/all_by_all_estimated_proteins_fkerg.pdf', 
       height = 9, width = 9)

# plotting estimated correlation matrix -----------------------------------
correlation_matrix_estimate_me <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(R[n,m]) %>% 
  median_qi() %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(n = number_val), 
             by = 'n') %>% 
  inner_join(sample_number_name_df_me %>% 
               dplyr::rename(m = number_val), 
             by = 'm')

heatmap_lower <- correlation_matrix_estimate_me %>% 
  ggplot(aes(x = protein_name_nice.x, y = protein_name_nice.y)) + 
  geom_raster(aes(fill = .lower)) + 
  xlab('') +
  ylab('') +
  ggtitle('A. 2.5th Percentile') +
  theme_bw() +
  scale_fill_gradient2(name = "Correlation")

heatmap_middle <- correlation_matrix_estimate_me %>% 
  ggplot(aes(x = protein_name_nice.x, y = protein_name_nice.y)) + 
  geom_raster(aes(fill = R)) + 
  xlab('') +
  ylab('') +
  ggtitle('B. 50th Percentile') +
  theme_bw() +
  scale_fill_gradient2(name = "Correlation")

heatmap_upper <- correlation_matrix_estimate_me %>% 
  ggplot(aes(x = protein_name_nice.x, y = protein_name_nice.y)) + 
  geom_raster(aes(fill = .upper)) + 
  xlab('') +
  ylab('') +
  ggtitle('C. 97.5th Percentile') +
  theme_bw() +
  scale_fill_gradient2(name = "Correlation")

heatmap_out_aggregated <- ggarrange(heatmap_lower, 
                                    heatmap_middle, 
                                    heatmap_upper, 
                                    nrow = 1, 
                                    common.legend = TRUE)

ggsave(heatmap_out_aggregated, filename = 'figures/heatmap_out_aggregated.svg', 
       width = 14, height = 4.5)


# generated quantities block ----------------------------------------------

## summary statistics for paper
estimate_pars_me_cov_10_1 %>% 
  spread_draws(mnfe_umol_per_mol_c) %>% 
  ## this needs to be truncated so filtering out values below zero
  dplyr::filter(mnfe_umol_per_mol_c > 0) %>% 
  median_qi()

## summary statistics for paper
estimate_pars_me_cov_10_1 %>% 
  spread_draws(cuzn_umol_per_mol_c) %>% 
  ## this needs to be truncated so filtering out values below zero
  dplyr::filter(cuzn_umol_per_mol_c > 0) %>% 
  median_qi()

## Figure of Posterior prediction
mnsod_post_10 <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(mnfe_umol_per_mol_c) %>% 
  ggplot(aes(x = mnfe_umol_per_mol_c)) +
  geom_histogram(fill = 'firebrick4', colour = 'firebrick4') +
  theme_bw() +
  xlim(0, 3) +
  ggtitle('C.') +
  ylab('Posterior Probability\nSample Count') +
  theme(plot.caption = element_text(hjust = -1, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(MnSOD~Contribution~to~Mn:C~(mu*mol:mol)));mnsod_post_10

cuznsod_post_10 <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(cuzn_umol_per_mol_c) %>% 
  ggplot(aes(x = cuzn_umol_per_mol_c)) +
  geom_histogram(fill = 'cadetblue4', colour = 'cadetblue4') +
  theme_bw() +
  ggtitle('D.') +
  xlim(0, 1) +
  ylab('') + 
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(CuZnSOD~Contribution~to~Cu*","~Zn:C~(mu*mol:mol)));cuznsod_post_10

mnsod_post_1 <- estimate_pars_me_cov_1_1 %>% 
  spread_draws(mnfe_umol_per_mol_c) %>% 
  ggplot(aes(x = mnfe_umol_per_mol_c)) +
  geom_histogram(fill = 'firebrick4', colour = 'firebrick4') +
  theme_bw() +
  xlim(0, 3) +
  ggtitle('H.') +
  ylab('Posterior Probability\nSample Count') +
  theme(plot.caption = element_text(hjust = -1, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(MnSOD~Contribution~to~Mn:C~(mu*mol:mol)));mnsod_post_1

cuznsod_post_1 <- estimate_pars_me_cov_1_1 %>% 
  spread_draws(cuzn_umol_per_mol_c) %>% 
  ggplot(aes(x = cuzn_umol_per_mol_c)) +
  geom_histogram(fill = 'cadetblue4', colour = 'cadetblue4') +
  theme_bw() +
  ggtitle('I.') +
  xlim(0, 1) +
  ylab('') + 
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(CuZnSOD~Contribution~to~Cu*","~Zn:C~(mu*mol:mol)));cuznsod_post_1

mnsod_post_2 <- estimate_pars_me_cov_2_1 %>% 
  spread_draws(mnfe_umol_per_mol_c) %>% 
  ggplot(aes(x = mnfe_umol_per_mol_c)) +
  geom_histogram(fill = 'firebrick4', colour = 'firebrick4') +
  theme_bw() +
  xlim(0, 3) +
  ggtitle('E.') +
  ylab('Posterior Probability\nSample Count') +
  theme(plot.caption = element_text(hjust = -1, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(MnSOD~Contribution~to~Mn:C~(mu*mol:mol)));mnsod_post_2

cuznsod_post_2 <- estimate_pars_me_cov_2_1 %>% 
  spread_draws(cuzn_umol_per_mol_c) %>% 
  ggplot(aes(x = cuzn_umol_per_mol_c)) +
  geom_histogram(fill = 'cadetblue4', colour = 'cadetblue4') +
  theme_bw() +
  ggtitle('F.') +
  xlim(0, 1) +
  ylab('') + 
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  xlab(expression(CuZnSOD~Contribution~to~Cu*","~Zn:C~(mu*mol:mol)));cuznsod_post_2

estimate_pars_me_cov_10_1 %>% 
  spread_draws(percentage_metallated) %>% 
  ggplot(aes(percentage_metallated)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw()

estimate_pars_me_cov_1_1 %>% 
  spread_draws(percentage_metallated) %>% 
  ggplot(aes(percentage_metallated)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw()

estimate_pars_me_cov_2_1 %>% 
  spread_draws(percentage_metallated) %>% 
  ggplot(aes(percentage_metallated)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw()

prior_21_p <- data.frame(sim_out = rbeta(n = 1000000, shape1 = 2, shape2 = 1)) %>% 
  ggplot(aes(sim_out)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw() +
  ylab('Probability Density') +
  xlab(expression(zeta~(Proportion~Protein~Metallated)))

prior_11_p <- data.frame(sim_out = rbeta(n = 1000000, shape1 = 1, shape2 = 1)) %>% 
  ggplot(aes(sim_out)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw() +
  ylab('Probability Density') +
  xlab(expression(zeta~(Proportion~Protein~Metallated)))

prior_101_p <- data.frame(sim_out = rbeta(n = 1000000, shape1 = 10, shape2 = 1)) %>% 
  ggplot(aes(sim_out)) +
  geom_density() +
  xlim(0, 1) +
  theme_bw() +
  ylab('Probability Density') +
  xlab(expression(zeta~(Proportion~Protein~Metallated)))

layer1 <- ggarrange(mnsod_post_10 + ggtitle('B.'), cuznsod_post_10 + ggtitle('C.'), prior_101_p + ggtitle('D.'), nrow = 1, align = 'hv')
layer1_main <- ggarrange(prior_101_p + ggtitle('B.'), mnsod_post_10, cuznsod_post_10, nrow = 1, align = 'hv')
layer3 <- ggarrange(mnsod_post_2, cuznsod_post_2, prior_21_p + ggtitle('G.'), nrow = 1, align = 'hv')
layer2 <- ggarrange(mnsod_post_1, cuznsod_post_1, prior_11_p + ggtitle('J.'), nrow = 1, align = 'hv')

stoichiometry_out <- ggarrange(ggplot() + theme_void() + ggtitle('A.'), 
                               layer1, layer3, layer2, heights = c(0.1, 0.3, 0.3, 0.3), nrow = 4)

stoichiometry_out_mostly_metal <- ggarrange(ggplot() + theme_void() + ggtitle('A.'), 
                                            layer1_main, heights = c(0.25, 0.75), nrow = 2)

ggsave(stoichiometry_out, filename = 'figures/stoichiometry_out_supp.svg', 
       height = 10, width = 14)
ggsave(stoichiometry_out_mostly_metal, filename = 'figures/stoichiometry_out_most_metal.svg', 
       height = 10*0.4, width = 14)

## Mean value
plot_of_estimated_mean_prot_per_carbon <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(ug_prot_per_ug_carbon_frag) %>% 
  ggplot(aes(ug_prot_per_ug_carbon_frag)) + 
  geom_histogram() +
  ylab('Posterior Probability\nSample Count') +
  xlab(expression(mu*g~Protein~mu*g~Carbon^-1)) +
  theme_bw();plot_of_estimated_prot_per_carbon

## SD
plot_of_estimated_sd_prot_per_carbon <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(sigma_protein_per_carbon_frag) %>% 
  ggplot(aes(sigma_protein_per_carbon_frag)) + 
  geom_histogram() +
  ylab('Posterior Probability\nSample Count') +
  xlab(expression(mu*g~Protein~mu*g~Carbon^-1)) +
  theme_bw();plot_of_estimated_sd_prot_per_carbon

posterior_prediction_for_protein_to_carbon <- estimate_pars_me_cov_10_1 %>% 
  spread_draws(ug_prot_per_ug_carbon_post_pred) %>% 
  ggplot(aes(ug_prot_per_ug_carbon_post_pred)) + 
  geom_histogram() +
  xlim(0, 2) +
  ylab('Posterior Probability Predictive Distribution\nSample Count') +
  xlab(expression(mu*g~Protein~mu*g~Carbon^-1)) +
  theme_bw();posterior_prediction_for_protein_to_carbon

## plotting the prior
rnorm(n = 100000, 
      mean = mean(liefer_data$protein_per_carbon), 
      sd = sd(liefer_data$protein_per_carbon)) %>% 
  as.data.frame() %>% 
  ggplot(aes(.)) +
  xlim(0, 2) +
  geom_histogram() +
  ylab('Posterior Probability Predictive Distribution\nSample Count') +
  theme_bw()

ggsave(posterior_prediction_for_protein_to_carbon, 
       filename = 'figures/prot_per_carbon_posterior.svg', width = 6.4, height = 5)

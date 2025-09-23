
## Making figures of basic biogeography.

library(reshape2)

all_avail_csvs_normed_cyl_kerg_abundance_w_peps_read_in <- read.csv('data/intermediate_data/all_avail_csvs_normed_cyl_kerg_abundance_w_peps.csv')
total_taxonomic_biomass_no_unambiguous_renorm_read_in <- read.csv('data/intermediate_data/total_taxonomic_biomass_no_unambiguous_renorm.csv')

## read in protein concentration data
prot_conc_data <- read.csv("data/protein_concentration_per_sample.csv")

df_mapping <- data.frame(assigned_tax_coarse = c(rep('Haptophyte', 2), rep('Diatom', 3)), 
                         assigned_tax = c('Other Haptophyta', 'Phaeocystis', 'Diatom', 'Fragilariopsis', 'Pseudo-nitzschia')) 

hapt_diatom_df_form <- total_taxonomic_biomass_no_unambiguous_renorm_read_in %>% 
  dplyr::filter(assigned_tax == 'Other Haptophyta' | assigned_tax == "Phaeocystis" | assigned_tax == 'Diatom' | assigned_tax == 'Fragilariopsis' | assigned_tax == 'Pseudo-nitzschia') %>%
  # dplyr::filter(assigned_tax == "Phaeocystis") %>% 
  inner_join(df_mapping, by = 'assigned_tax') %>% 
  group_by(sample_id, assigned_tax, assigned_tax_coarse) %>% 
  summarize(mean_summed_prot_biomass_no_unk_renorm = mean(summed_prot_biomass_no_unk_renorm)) %>% 
  group_by(sample_id, 
           assigned_tax_coarse) %>%
  summarize(sum_coarse_prot_biomass = sum(mean_summed_prot_biomass_no_unk_renorm)) %>% 
  dcast(sample_id ~ assigned_tax_coarse) %>% 
  ## join it with the protein concentration data
  inner_join(prot_conc_data %>% 
               dplyr::rename(sample_id = sample_id_nlz), by = 'sample_id') %>% 
  dplyr::mutate(diatom_prot_conc = Diatom*protein_conc_ug_l,
                hapto_prot_conc = Haptophyte*protein_conc_ug_l)

hapt_diatom_plot <- hapt_diatom_df_form %>% 
  ggplot(aes(x = diatom_prot_conc, y = hapto_prot_conc)) +
  geom_point(size = 4, aes(fill = diatom_prot_conc), pch = 21) +
  scale_fill_gradient(high = "darkgreen", low = "white", 
                      limits = c(0.1, 50), 
                      name = expression("Diatom-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")")) +
  theme_bw() +
  # theme(legend.position = 'none') +
  theme(legend.position = c(0.7, 0.7), 
        legend.direction = "horizontal",
        legend.box.background = element_rect(color = "grey40", size = 0.6), # Add black border
  ) +
  xlab(expression("Diatom-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")")) +
  ylab(expression("Haptophyte-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")"));hapt_diatom_plot

hapt_diatom_plot_legend_only <- hapt_diatom_df_form %>% 
  ggplot(aes(x = diatom_prot_conc, y = hapto_prot_conc)) +
  geom_point(size = 4, aes(fill = diatom_prot_conc), pch = 21) +
  scale_fill_gradient(high = "darkgreen", low = "white", 
                      limits = c(0.1, 50), 
                      name = expression("Diatom-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")")) +
  theme_classic() +
  # theme(legend.position = 'none') +
  theme(legend.position = c(0.5, 0.5), 
        legend.direction = "horizontal",
        legend.box.background = element_rect(color = "grey40", size = 0.6), # Add black border
  ) +
  xlim(0, 100) +
  ylim(0, 200) +
  xlab(expression("Diatom-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")")) +
  ylab(expression("Haptophyte-associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")"));hapt_diatom_plot_legend_only

## plot F. kerguelensis
kerg_cyl_plot <- all_avail_csvs_normed_cyl_kerg_abundance_w_peps_read_in %>% 
  group_by(sample_id) %>% 
  summarize(kerg_mean = mean(fraction_kerg),
            cyl_mean = 1 - mean(fraction_kerg)) %>% 
  inner_join(hapt_diatom_df_form, by = 'sample_id') %>% 
  dplyr::mutate(kerg_prot_conc = Diatom*kerg_mean*protein_conc_ug_l,
                cyl_prot_conc = Diatom*cyl_mean*protein_conc_ug_l) %>% 
  ggplot(aes(x = kerg_prot_conc, y = cyl_prot_conc)) +
  geom_point(size = 4, aes(fill = diatom_prot_conc), pch = 21) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 50) +
  ylim(0, 50) +
  theme_bw() +
  coord_equal() +
  scale_fill_gradient(high = "darkgreen", low = "white", limits = c(0.1, 50)) +
  theme(legend.position = 'none') +
  xlab(expression(italic("F. kerguelensis")~"associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")")) +
  ylab(expression(italic("F. cylindrus")~"associated Protein Concentration (" * mu * "g" ~ L^{-1} * ")"));kerg_cyl_plot

## make a dataframe that maps the station id to the sample id
sample_id_to_station_map <- env_filter_data %>% 
  dplyr::select(sample_id, station) %>% 
  unique()

mean_taxa_by_station <- hapt_diatom_df_form %>% 
  dplyr::mutate(sample_id = as.character(sample_id)) %>% 
  inner_join(sample_id_to_station_map, by = 'sample_id') %>% 
  group_by(station) %>% 
  summarize(mean_diatom_by_station = mean(diatom_prot_conc),
            mean_haptophyte_by_station = mean(hapto_prot_conc))

source('scripts/sample_map.R')

composite_of_map_and_profiles <- ggarrange(ggarrange(map_of_sampling + ggtitle('A.'), 
                    ggplot() + theme_void() + ggtitle('         B.'), 
                    widths = c(0.5, 0.5)),
          ggarrange(hapt_diatom_plot + ggtitle('C.'), 
                    profiles_proteins_together + ggtitle('D.')), nrow = 2);composite_of_map_and_profiles
                                           
ggsave(composite_of_map_and_profiles, filename = 'figures/composite_of_map_and_figures.pdf', width = 8, height = 8)
ggsave(kerg_cyl_plot, filename = 'figures/kerguelensis_vs_cylindrus_concentrations.pdf', height = 7, width = 7)
ggsave(hapt_diatom_plot_legend_only, filename = "figures/diatom_colour_bar.pdf")

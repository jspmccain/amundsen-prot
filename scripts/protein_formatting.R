# matching metal concentrations with proteins

library(readxl)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)

# reading in data and formatting.

## This file has station, casts, trace metal concentrations, and other environmental variables.
amund_data <- read_excel('data/Data selection_ASP.xlsx',
                         sheet = 1)

# changing the column name format
names(amund_data) <- gsub(pattern = " ", replacement = "_", 
                          x = tolower(names(amund_data)),
                          fixed = TRUE)

amund_data_full <- amund_data[grepl(pattern = "ANA", x = amund_data$cruise),]

### getting protein filter information. This file has station, depths, bottle numbers, and depths recorded for protein samples. These are *target* depths that were manually recorded at sea.
prot_filters <- read_excel("data/protein_sample_sheets.xlsx",
                           sheet = 1)

## append the actual par values
light_data_from_ctd <- read_excel(path = 'data/merged bottle files all data in one sheet_RobMiddag_2023.xlsx', 
                                  sheet = 1) %>% 
  dplyr::rename(PAR_adjusted = PAR,
                bottle = Bottle)

## extract the station number from the station_bottle string
station_as_numeric <- substr(light_data_from_ctd$Station_cast, 1, 2) %>% as.numeric()
light_data_from_ctd$station <- station_as_numeric

## appending the depths fired and the adjusted par levels
prot_filters_appended <- prot_filters %>% 
  dplyr::rename(target_depth = depth,
                bottle = niskin_number) %>% 
  inner_join(light_data_from_ctd, by = c('station', 'bottle')) %>% 
  ## getting the PAR and temperature for the same bottles as the protein bottles
  dplyr::rename(par_from_protein_bottles = PAR_adjusted,
                temp_from_protein_bottles = Temperature)
 
## now make a dataframe that matches depths
make_depth_matched_df <- function(){
  
  ## make an empty dataframe with all the correct column names
  empty_df <- amund_data_full[0, ]
  
  ## loop through all the protein filters
  for(row_i in 1:nrow(prot_filters_appended)){
    
    ## get the station name from the protein filter dataframe, and the depth
    station_i <- prot_filters_appended[row_i, ]$station
    depth_i <- prot_filters_appended[row_i, ]$Depth
    
    ## filter metal data based on station
    station_metal_data <- amund_data_full %>% 
      dplyr::filter(station == station_i)
    
    ## after filtering the station, find the closest depth.
    closest_corresponding_depth <- station_metal_data[which.min(abs(depth_i - station_metal_data$`depth_[m]`)),]
    
    # Some of the protein filters do not have metal data
    if(nrow(closest_corresponding_depth) == 0){
      
      ## this is just to format the columns
      closest_corresponding_depth <- station_metal_data[1,]
      
      ## this makes all the column values NA
      closest_corresponding_depth[1,] <- NA
    }
    
    empty_df <- rbind(empty_df, 
                      closest_corresponding_depth)
    
  }
  return(empty_df)
}

depth_matched_df <- make_depth_matched_df()

## after this depth matched dataframe has been made, we do a column bind.
merged_metal_prot <- cbind(prot_filters_appended, depth_matched_df %>% rename(depth_m_tm_bottle = `depth_[m]`,
                                                             cast_tm_bottle = cast,
                                                             station_tm_bottle = station,
                                                             salinity_tm_bottle = salinity))

write.csv(merged_metal_prot %>% 
            dplyr::rename(temperature_c = `temperature_[øc]`,
                          sbeox = sbeox0ps,
                          oxygen_mol_per_kg_check = `oxygen_mol_[æmol/kg]`,
                          oxygen_molar = `oxygen_molar_[æm]`), 
          file = "data/protein_sample_sheets_with_metal_data.csv")

## This spreadsheet has a new column and is modified from the above one. It now has a manual selection, so that file name is appended with this string.
## The column "depth_m" corresponds with "Depth" in the sheet "protein_sample_sheets_with_metal_data.csv". This is the depth at which the protein bottle fired, and is used later to plot.
man_sel <- read.csv("data/protein_sample_sheets_with_metal_data_manual_selection.csv", check.names = FALSE)

#making a new column with just the number
man_sel$sample_id <- str_sub(man_sel$sample_label, start = -3)

### joining protein amounts to the manual selection sheet
protein_amounts <- read.csv("data/Amundsen2018_2018_BCA_all_pages.csv")

## reformatting the sample _ID column in include preceding zeros
include_preceding_zeros <- function(protein_amounts_df){
  protein_amounts_df$sample_id <- as.character(protein_amounts_df$sample_id) 
  for(i in 1:nrow(protein_amounts_df)){
    if(nchar(protein_amounts_df[i, ]$sample_id) == 2){
      protein_amounts_df[i, ]$sample_id <- paste0('0', protein_amounts_df[i, ]$sample_id)
    }
  }
  return(protein_amounts_df)
}

#formatted protein amounts so that the sample_id matches the metals above.
protein_amounts_form <- include_preceding_zeros(protein_amounts_df = protein_amounts)

# joining the protein amounts with the sample sheet + metal concentrations
man_sel_with_protein <- man_sel %>% 
  inner_join(protein_amounts_form %>% 
                         dplyr::select(sample_id, total_prot_amount), 
                                                    by = 'sample_id')

sample_chose <- man_sel_with_protein %>% dplyr::filter(sample_choices == 'Y')

mn_selections_p <- man_sel_with_protein %>%
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(x = depth_m_tm_bottle, y = `dmn_[nm]`)) +
  geom_point() +
  ylab('dMn (nM)') +
  coord_flip() +
  xlab('Depth (m)') +
  # ylim(0, 1) +
  theme_bw() +
  facet_wrap(~station) +
  geom_line() +
  scale_x_reverse() +
  geom_vline(data = sample_chose, aes(xintercept = depth_m));mn_selections_p

fe_selections_p <- man_sel_with_protein %>%
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(x = depth_m_tm_bottle, y = `dfe_[nm]`)) +
  geom_point() +
  ylab('dFe (nM)') +
  xlab('Depth (m)') +
  coord_flip() +
  ylim(0, 1) +
  theme_bw() +
  facet_wrap(~station) +
  geom_line() +
  scale_x_reverse() +
  geom_vline(data = sample_chose, aes(xintercept = depth_m));fe_selections_p

light_selections_p <- merged_metal_prot %>%
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(x = Depth, y = par_from_protein_bottles)) +
  ylab(expression("Photosynthetically Active Radiation (" * mu * "mol photons " * m^-2 * " sec"^-1 * ")")) +
  xlab('Depth (m)') +
  geom_point() +
  coord_flip() +
  # ylim(0, 1) +
  # xlim(0, 50) +
  theme_bw() +
  facet_wrap(~station) +
  geom_line() +
  scale_x_reverse() +
  geom_vline(data = sample_chose, aes(xintercept = depth_m));light_selections_p

temp_selections_p <- man_sel_with_protein %>%
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(x = depth_m, y = temperature)) +
  ylab('Temperature (C)') +
  xlab('Depth (m)') +
  geom_point() +
  coord_flip() +
  # ylim(0, 1) +
  theme_bw() +
  facet_wrap(~station) +
  geom_line() +
  scale_x_reverse();temp_selections_p

ggsave(light_selections_p, filename = 'figures/light_selections.pdf', width = 7, height = 6)  
ggsave(fe_selections_p, filename = 'figures/fe_selections.pdf', width = 7, height = 6) 
ggsave(mn_selections_p, filename = 'figures/mn_selections.pdf', width = 7, height = 6)
ggsave(temp_selections_p, filename = 'figures/temp_profiles.pdf', width = 7, height = 6)

man_sel_with_protein %>% 
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(y = total_prot_amount/volume_filtered, x = target_depth)) +
  ylab('Protein Amounts (ug / L)') +
  geom_point() +
  xlab('Depth (m)') +
  coord_flip() +
  # ylim(0, 1) +
  theme_bw() +
  facet_wrap(~station) +
  geom_line() +
  scale_x_reverse()

man_sel_with_protein %>% 
  dplyr::filter(station != 1, station != 4, station != 10, station != 11) %>% 
  ggplot(aes(y = total_prot_amount/volume_filtered , x = `pon_concentration_mg/ml`)) +
  geom_point(size = 3) +
  # scale_x_log10() +
  theme_bw() +
  ylab('Protein Concentration (ug/L)') +
  xlab('PON Concentration (mg/mL)') +
  geom_label(aes(label = station), nudge_y = 3)

### I want to format a sample_id column from the protein dataframe so that I can bind it with the metaproteomic fractions. This needs all the hanging zeros to be removed.
remove_leading_zero <- function(strings) {
  # Use gsub to check if the string starts with a zero and remove it
  result <- gsub("^0", "", strings)
  return(result)
}

man_sel_with_protein$sample_id_nlz <- remove_leading_zero(man_sel_with_protein$sample_id)

# Make a dataframe that has the protein concentration so that I can write it out and use it elsewhere.
protein_concentration_df <- man_sel_with_protein %>% 
  dplyr::mutate(protein_conc_ug_l = total_prot_amount/volume_filtered) %>% 
  dplyr::select(sample_id_nlz, protein_conc_ug_l)

protein_concentration_df_st_targ <- man_sel_with_protein %>% 
  dplyr::mutate(protein_conc_ug_l = total_prot_amount/volume_filtered) %>% 
  dplyr::select(sample_id_nlz, protein_conc_ug_l, station, target_depth, total_prot_amount, volume_filtered)

write.csv(x = protein_concentration_df, file = 'data/protein_concentration_per_sample.csv', row.names = FALSE)
write.csv(x = protein_concentration_df_st_targ, file = 'data/protein_concentration_per_sample_station_target_depth.csv', row.names = FALSE)


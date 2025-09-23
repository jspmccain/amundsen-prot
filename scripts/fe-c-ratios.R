# plotting particulate fe to carbon ratios

# matching metal concentrations with proteins

library(readxl)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)

# reading in data and formatting
amund_data <- read_excel('data/Data selection_ASP.xlsx',
                         sheet = 1)

amund_data_transformed <- amund_data %>% 
  rowwise() %>% 
  dplyr::mutate(pfe_labile_mol_l = `PFe [nM] labile`*1e-3, ## converts to umol / L (Converting from nMol to uMol = 1e6/1e9)
                c_mol_per_l = `POC concentration mg/ml`*(1/12), ## converts to mol C per L
                pfe_umol_per_mol_c = pfe_labile_mol_l/c_mol_per_l,
                pmn_labile_mol_l = `PMn [nM] labile`*1e-3,
                pmn_umol_per_mol_c = pmn_labile_mol_l/c_mol_per_l)

depth_profiles_of_p_mn <- amund_data_transformed %>% 
  # dplyr::filter(Station == 50) %>% dplyr::select('PMn [nM] labile')
  dplyr::filter(Station %in% c(24, 52, 55, 57)) %>% 
  dplyr::filter(`Depth [m]` < 40) %>%
  ggplot(aes(x = pmn_umol_per_mol_c, y = `Depth [m]`)) +
  geom_point() +
  # scale_x_log10() +
  scale_y_reverse() + 
  facet_wrap(~Station) +
  # coord_flip() +
  xlab(expression("Particulate Mn " * mu * "mol / Mol C")) +
  theme_bw() +
  ggtitle(expression("Mn:C (" * mu * "mol:mol), Depth < 40m"));depth_profiles_of_p_mn

ggsave(depth_profiles_of_p_mn, filename = 'figures/depth_profiles_of_p_mn.pdf', height = 4, width = 5)


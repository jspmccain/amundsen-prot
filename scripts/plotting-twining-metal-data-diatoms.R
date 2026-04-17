### looking at Twining metal quota data

library(dplyr)
library(magrittr)
library(readxl)

twining_data <- read.csv("data/MnC_SXRF_data compilation_22Jan2021_forScott.csv")

twining_data_diatoms <- twining_data %>% 
  dplyr::filter(CellType2 == "diatom")

twining_data_diatoms_c_mn <- complete.cases(twining_data_diatoms[,c("cellC", "cellMn")])

median_value_mn_c <- median(c(1000000*twining_data_diatoms$cellMn/twining_data_diatoms$cellC), na.rm = TRUE)

mn_c_ratio_twining_plot <- twining_data_diatoms[twining_data_diatoms_c_mn, ] %>% 
  ggplot(aes(x = 1000000*cellMn/cellC)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  ggtitle('A. Single cell diatom Mn:C') +
  xlab(expression(mu*"mol Mn : mol C")) +
  ylab('Count') +
  annotate(geom = "label", x = 0.1, y = 35, label = "n = 433") +
  geom_vline(xintercept = median_value_mn_c,
             lty = 2);mn_c_ratio_twining_plot

ggsave(plot = mn_c_ratio_twining_plot, filename = 'figures/mn_c_ratio_single.pdf')


# plotting only southern ocean --------------------------------------------

twining_data_diatoms_southern_ocean <- twining_data_diatoms %>% 
  dplyr::filter(Lat_N < -30)

twining_data_diatoms_southern_ocean_cc_mn_c <- complete.cases(twining_data_diatoms_southern_ocean[,c("cellC", "cellMn")])

median_value_mn_c_southern_ocean <- median(c(1000000*twining_data_diatoms_southern_ocean$cellMn/twining_data_diatoms_southern_ocean$cellC), na.rm = TRUE)

mn_c_ratio_twining_plot_southern_ocean <- twining_data_diatoms_southern_ocean[twining_data_diatoms_southern_ocean_cc_mn_c, ] %>% 
  ggplot(aes(x = 1000000*cellMn/cellC)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  ggtitle('B. Single cell diatom Mn:C (Southern Ocean only)') +
  xlab(expression(mu*"mol Mn : mol C")) +
  ylab('Count') +
  annotate(geom = "label", x = 0.3, y = 5, label = "n = 63") +
  geom_vline(xintercept = median_value_mn_c,
             lty = 2);mn_c_ratio_twining_plot_southern_ocean



# adding them together ----------------------------------------------------

ggsave(plot = ggarrange(mn_c_ratio_twining_plot, mn_c_ratio_twining_plot_southern_ocean, nrow = 2),
       filename = 'figures/mn_c_ratio_split.pdf')

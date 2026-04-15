# setwd('~/My Drive (jspmccain@gmail.com)/Projects/amundsen-prot/')

## This script goes through the Fragilariopsis cylindrus proteomic data and makes figures.

library(magrittr)
library(dplyr)
library(GGally)
library(reshape2)
library(ggplot2)
library(readr)
library(ggpubr)
library(ggExtra)

## A dataframe that maps the replicate number to the treatment.
df_treat <- data.frame(replicate = c("01", "02", "03","04", "05", "06", "07", "08", "09", "10", "11", "12"), 
                       treat = c(rep(c('Low Fe', 'Medium Fe', 'High Fe', 'Low Mn, High Fe'), each = 3)))
df_treat$replicate <- as.factor(df_treat$replicate)

## Order the factors to help the plots look nicer later on.
df_treat$treat <- factor(df_treat$treat, levels = c('Low Fe', 'Medium Fe', 'High Fe', 'Low Mn, High Fe'))

## Read in the proteomic data.
prot_data <- read_tsv('data/dia_frag_iq_02.tsv') %>% 
  dplyr::rename(protein_data = Protein.Group,
                protein_list = `All Mapped Proteins`)

## Dimensions and characteristics of the proteomic data.
dim(prot_data)
head(prot_data)

## Looking that the names are set correctly.
head(prot_data)

## I want to do some targeted analysis of various proteins. Format the dataframe to be easy to work with in ggplot2.
prot_data_long_all <- prot_data %>% 
  dplyr::select(-n_fragments, -n_peptides) %>% 
  # dplyr::select(-n_fragments, -n_peptides, -all_mapped_prots) %>% 
  melt(variable.name = 'replicate') %>% 
  inner_join(df_treat, by = 'replicate')

## Looking at the proteins that have the word photosystem
prot_data_subset <- prot_data %>% 
  dplyr::filter(grepl(pattern = 'photosystem', x = protein_data))

# write.csv(prot_data_subset, 
#           file = 'data/frag_cyl_culture_data/photosystem_protein_data_frag_cyl_subset.csv',
#           row.names = FALSE)
# 
# ## manually grouped, and now reading in again
# prot_data_subset_ps <- read.csv('data/frag_cyl_culture_data/photosystem_protein_data_frag_cyl_subset_manual_groupings.csv')

# prot_data_subset_ps_long <- prot_data_subset_ps %>% 
#   dplyr::select(-n_fragments, -n_peptides) %>% 
#   melt(variable.name = 'replicate') %>% 
#   inner_join(df_treat, by = 'replicate')

## This is for making the axes pretty.
scaleFUN <- function(x) sprintf("%.1f", x)

## plotting the photosystem action cadetblue4
## Making a figure of every single protein in the photosystems
make_protein_plot <- function(protein_data_name, df_input = prot_data_subset_ps_long, 
                              custom_title = NA, 
                              colour_choice = 'firebrick4',
                              custom_axis_label = 'Protein Abundance'){
  ## filter filter for a specific protein naem
  # df_input <- prot_data_subset_ps_long
  # protein_data_name <- prot_data_subset_ps$protein_data[1]
  # 
  
  ## filter the dataframe for your protein of interest
  df_subset <- df_input %>% 
    dplyr::filter(protein_data == protein_data_name)
  
  ## make the figure
  output_fig <- df_subset %>% 
    ggplot(aes(x = treat, y = value)) +
    geom_hline(yintercept = df_subset$value %>% max(na.rm = TRUE),
               lty = 2, colour = 'grey70') + 
    geom_hline(yintercept = df_subset$value %>% min(na.rm = TRUE),
               lty = 2, colour = 'grey70') +
    geom_point(size = 3, fill = colour_choice, pch = 21) +
    theme_bw() +
    ggtitle(protein_data_name) +
    xlab('Treatment') +
    # theme(plot.margin = unit(c(0,1.5,0,0), "cm")) +
    ylab(custom_axis_label) +
    scale_y_continuous(labels = scaleFUN) +
    theme(plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  
  if(!is.na(custom_title)){
    output_fig <- output_fig + ggtitle(custom_title)
  }
  
  output_fig_hist <- prot_data_long_all %>% 
    dplyr::filter(replicate == "01") %>% 
    ggplot(aes(value)) +
    geom_density() +
    coord_flip() +
    theme_bw() +
    ylab('Density') +
    # xlab('') +
    labs(x = NULL) +
    theme(
      plot.margin = margin(t = 0,  # Top margin
                                         r = 0,  # Right margin
                                         b = 0,  # Bottom margin
                                         l = 0),  # Reduce margin
      axis.title.y = element_blank()
    ) +
    geom_vline(xintercept = df_subset$value %>% max(na.rm = TRUE),
               lty = 2) +
    geom_vline(xintercept = df_subset$value %>% min(na.rm = TRUE), 
               lty = 2) + 
    annotate("rect", fill = colour_choice, alpha = 0.5, 
             xmin = df_subset$value %>% min(na.rm = TRUE), 
             xmax = df_subset$value %>% max(na.rm = TRUE),
             ymin = -Inf, ymax = Inf)
  
  arranged_out <- ggarrange(output_fig, output_fig_hist, widths = c(0.75, 0.25), 
                            align = 'hv')
  
  return(arranged_out)
}

make_several_protein_plots <- function(list_of_protein_names, df_input){
  ## loop through the list of protein names
  for(i in 1:length(list_of_protein_names)){
    p_out <- make_protein_plot(protein_data_name = list_of_protein_names[i], 
                      df_input = df_input)
    
    ## if the protein name has a slash in it, remove
    protein_name_i_no_slash <- gsub("/", "", list_of_protein_names[i])
    
    
    ggsave(plot = p_out, filename = paste0('figures/frag_cyl_prot/', protein_name_i_no_slash, '.pdf'), 
           height = 5, width = 7)
  }
}

# make_several_protein_plots(list_of_protein_names = unique(prot_data_subset_ps_long$protein_data), df_input = prot_data_subset_ps_long)

## now I want to examine superoxide dismutase changes
prot_data_superoxide <- prot_data %>% 
  dplyr::filter(grepl(pattern = 'superoxide', x = protein_data)) %>% 
  dplyr::select(-n_fragments, -n_peptides) %>% 
  melt(variable.name = 'replicate') %>% 
  inner_join(df_treat, by = 'replicate')

# Trying to synthesize SOD plots ------------------------------------------

nisod <- make_protein_plot(protein_data_name = unique(prot_data_superoxide$protein_data)[1], 
                  df_input = prot_data_superoxide, custom_title = 'D. NiSOD (OEU22085.1)',
                  colour_choice = 'lightgoldenrod', custom_axis_label = '')
cuznsod1 <- make_protein_plot(protein_data_name = unique(prot_data_superoxide$protein_data)[2], 
                  df_input = prot_data_superoxide, custom_title = 'C. CuZnSOD (OEU15804.1)', 
                  colour_choice = 'cadetblue4')
mnsod1 <- make_protein_plot(protein_data_name = unique(prot_data_superoxide$protein_data)[3], 
                  df_input = prot_data_superoxide, custom_title = 'A. MnSOD (OEU16864.1)')
mnsod2 <- make_protein_plot(protein_data_name = unique(prot_data_superoxide$protein_data)[4], 
                  df_input = prot_data_superoxide, custom_title = 'B. MnSOD (OEU17064.1)', custom_axis_label = '')
# cuznsod2 <- make_protein_plot(protein_data_name = unique(prot_data_superoxide$protein_data)[5], 
#                   df_input = prot_data_superoxide, custom_title = 'B. CuZnSOD (OEU16864.1)', 
#                   colour_choice = 'cadetblue4')

sod_regulation_figures <- ggarrange(mnsod1, mnsod2, cuznsod1, nisod)

ggsave(sod_regulation_figures, filename = 'figures/sod_regulation_figures.svg', width = 11, height = 7)



# Make a correlation plot for all the SODs --------------------------------

library(GGally)

protein_df_nice_names <- data.frame(protein_data = prot_data_superoxide$protein_data %>% unique(),
                                    protein_data_nn = c('NiSOD (OEU22085.1)', 
                                                        'CuZnSOD (OEU15804.1)', 
                                                        'MnSOD (OEU16864.1)', 
                                                        'MnSOD (OEU17064.1)', 
                                                        'CuZnSOD (OEU19597.1)'))

all_by_all_sod_scatterplot <- prot_data_superoxide %>%
  inner_join(protein_df_nice_names, by = 'protein_data') %>% 
  dplyr::select(-protein_data) %>% 
  dcast(replicate ~ protein_data_nn) %>% 
  dplyr::select(-replicate) %>% 
  ggpairs() +
  theme_bw() +
  xlab('Protein Abundance') +
  ylab('Protein Abundance')

ggsave(all_by_all_sod_scatterplot, filename = 'figures/all_by_all_sod_scatterplot.pdf',
       width = 13, height = 11)

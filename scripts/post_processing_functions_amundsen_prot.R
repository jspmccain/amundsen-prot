### openms post processing functions

parse_contigs <- function(protein_quant_out, row_i){
  
  ## input is the quant output from the FeatureFinderIdentification csv ProteinQUantifier (OpenMS)
  ## output is the list of ORFS associated with a given peptide in row_i
  
  contigs_prot <- protein_quant_out[row_i,]$protein %>% 
    as.character() %>% 
    strsplit(split = "/", fixed = TRUE) %>% 
    unlist()
  
  return(contigs_prot)
  
}

read_and_clean_file <- function(pep_file_name){
  # read in file with peptide mappings to contigs with peptide intensities
  prot_quant_output <- read.table(pep_file_name, 
                                  col.names = c('peptide', 
                                                'protein', 
                                                'n_proteins', 
                                                'charge', 
                                                'abundance'), 
                                  skip = 3)
  prot_quant_output_no_contam <- prot_quant_output[!grepl(pattern = "|", 
                                                          x = prot_quant_output$protein, 
                                                          fixed = TRUE), ]
  return(prot_quant_output_no_contam)
}

# functions for uploading data and manipulating files

get_multiple_files <- function(list_of_files){
  ## loops through read and clean file to then import all csv files from OpenMS output

  # read one  in so we don't have to format the empty_df
  empty_df <- data.frame(peptide = factor(),
                       protein = factor(),
                       n_proteins = numeric(),
                       charge = numeric(),
                       abundance = numeric(),
                       sample_id = numeric(),
                       file_string = character())
  
  # loop through all the rest of the files
  for(file_i in 1:length(list_of_files)){
    print(file_i)
    file_i_df <- read_and_clean_file(pep_file_name = list_of_files[file_i])
    # parsing the file name to retrieve the sample number
    file_i_df$sample_id <- rep(strsplit(list_of_files[file_i], "_")[[1]][5], nrow(file_i_df))
    file_i_df$file_string <- rep(list_of_files[file_i], 
                                 nrow(file_i_df))
    
    # append the file
    empty_df <- rbind(empty_df, file_i_df)
  }
  
  return(empty_df)

}

get_all_corresponding_proteins <- function(protein_quant_out){
  ## get all proteins corresponding to the peptides identified
  corresponding_proteins <- c()
  # loop through all the rows
  for(i in 1:nrow(protein_quant_out)){
    parse_contigs_i <- parse_contigs(protein_quant_out = protein_quant_out, 
                                     row_i = i)
    corresponding_proteins <- c(corresponding_proteins, 
                                parse_contigs_i)
  }
  
  return(corresponding_proteins)
  
}

get_tfg_peps <- function(tfg_ids, openms_quant_output){
  ## function to get all the peptides that correspond to tfg antioxidants
  
  empty_df <- data.frame(peptide = character(),
                         protein = character(),
                         n_proteins = numeric(),
                         charge = numeric(),
                         abundance = numeric())
  
  for(i in 1:length(tfg_ids)){
    
    sub_df_i <- openms_quant_output %>% filter(grepl(tfg_ids[i], 
                                                     protein))
    
    empty_df <- rbind(empty_df, sub_df_i)
  }
  
  return(empty_df)
}

normalize_out <- function(rep_file_list){
  
  print('Running normalization...')
  
  #########
  # Normalizes files across replicates
  #########
  
  # read in files with MCL assignments
  file_list <- list()
  for(file in 1:length(rep_file_list)){
    # file_list[[file]] <- read.csv(rep_file_list[file])
    file_list[[file]] <- read_and_clean_file(rep_file_list[file])
  }
  
  # make peptide list per file for determining common set of peptides across runs
  pep_list <- list()
  
  for(file in 1:length(file_list)){
    pep_list[[file]] <- file_list[[file]]$peptide
  }
  
  # find the common set of peptides across runs
  common_peps <- Reduce(intersect, pep_list)
  
  common_file_list <- list()
  
  # subset only common-across-runs peptides (to compare across samples), 
  # calculate db-dependent norm factor, apply norm factor
  for(file in 1:length(rep_file_list)){
    
    temp_df <- file_list[[file]]
    temp_df2 <- temp_df[temp_df$peptide %in% common_peps,]
    temp_df2$peptide <- temp_df2$peptide %>% as.character()
    
    # exclude contaminants from normalization factor calculation
    temp_df_no_contam <- temp_df[!grepl(pattern = "|", 
                                        x = temp_df$protein,
                                        fixed = TRUE), ]
    
    # Normalization factor calculated from all peptides identified (+ not common across MS runs)
    db_dep_norm_factor <- sum(temp_df_no_contam$abundance)
    temp_df2$db_norm_abund <- temp_df2$abundance/db_dep_norm_factor
    common_file_list[[file]] <- temp_df2
  }
  
  # making two matrices that have no normalization, and that has normalization  
  no_normal_abund_matrix <- matrix(nrow = length(common_peps), 
                                   ncol = length(rep_file_list))
  abund_matrix <- matrix(nrow = length(common_peps), 
                         ncol = length(rep_file_list))
  
  for(file in 1:length(rep_file_list)){
    temp_df <- common_file_list[[file]]
    no_normal_abund_matrix[, file] <- temp_df$abundance
    abund_matrix[ ,file] <- temp_df$db_norm_abund
  }
  
  # getting the column of names which will later be used.
  temp_file_no_abund <- dplyr::select(.data = common_file_list[[1]], 
                                      -c(db_norm_abund, abundance))
  
  # this does not change to the sum, it averages the abundance across
  abund_no_norm_avg <- rowMeans(no_normal_abund_matrix)
  db_abund_avg <- rowMeans(abund_matrix)
  
  consensus_file <- cbind(temp_file_no_abund, db_abund_avg, abund_no_norm_avg)
  
  return(consensus_file)
}

get_tax_assoc <- function(quant_out){
  ## function to assign taxonomy to peptides
  empty_tax_string <- c()
  # go through every peptide
  for(row_i in 1:nrow(all_anti_peps)){
    # if there is only one matched protein, take the taxonomic string from that protein
    if(all_anti_peps$n_proteins[row_i] == 1){
      # first convert it to a character      
      protein_i <- as.character(all_anti_peps[row_i, ]$protein)
      # if it was matched in the TARA set, then used the taxa string
      if(grepl(pattern = "TARA", x = protein_i)){
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        tax_string_i <- tara_tax[tara_tax$`Genome_Id final names` == protein_i_taxa, ]$tax_string
        tax_i <- tax_string_i
      } else {
        # take the taxa stringfrom the TFG annotation file
        tax_i <- annot_mcl[annot_mcl$orf_id == protein_i, ]$best_tax_string %>% as.character()
      }
      empty_tax_string <- c(empty_tax_string, tax_i)
      
    }
    # if it matched more than one protein, the plot thickens.
    if(all_anti_peps$n_proteins[row_i] > 1){
      # get the list of proteins
      protein_i <- parse_contigs(protein_quant_out = all_anti_peps, row_i = row_i)
      # this sum might look weird, it's because protein_i is a vector and the grepl statement takes a vector of TRUE/FALSE. The sum ofa 
      # vector of TRUE/FALSE is 1./0 respectively
      if(sum(grepl(pattern = "TARA", x = protein_i)) > 0){
        
        # taking the taxonomic string out of the protein_id
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        
        all_classes <- tara_tax[tara_tax$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$Best_taxonomy_CLASS
        
        if(length(unique(all_classes)) == 1){
          class_i <- unique(all_classes) %>% as.character()
          
          print(all_classes)
        } else {
          class_i <- 'non_unique_to_class'
        }
        
      } else {
        class_i <- paste(annot_mcl[annot_mcl$orf_id %in% protein_i, ]$best_tax_string %>% as.character(), collapse = ',')
      }
      
      empty_tax_string <- c(empty_tax_string, class_i)
    }
  }
  return(empty_tax_string)
}

classify_using_tfg <- function(tfg_annot, protein_i){
  ## small function for assigning taxonomic redundancy
  all_grpnorm_taxgrps <- tfg_annot[tfg_annot$orf_id %in% protein_i, ]$grpnorm_taxgrp %>% unique()
  if(length(all_grpnorm_taxgrps) == 1){
    class_i <- all_grpnorm_taxgrps 
  }
  else if(length(all_grpnorm_taxgrps) > 1){
    class_i <- "non-unique-grpnorm-taxgrp" 
  }
  else if(length(all_grpnorm_taxgrps) == 0){
    class_i <- "no-affiliated-grpnorm-taxgrp" 
  }
  else if(is.na(all_grpnorm_taxgrps)){
    class_i <- "no-affiliated-grpnorm-taxgrp" 
  }
  
  return(class_i)
}

classify_frag_using_tfg <- function(tfg_annot, protein_i, class_i){
  ## small function for assigning taxonomic redundancy
  if(class_i == 'Fragilariopsis'){
    all_tax_strings <- tfg_annot[tfg_annot$orf_id %in% protein_i, ]$best_tax_string %>% unique()
    if(length(all_tax_strings) == 1){
      tax_string_i <- all_tax_strings
    }
    ## make a new taxonomic string if they are not all unique
    else if(length(all_tax_strings) > 1){
      tax_string_i <- paste(all_tax_strings, collapse = '___')
    }
  } else if(class_i != 'Fragilariopsis'){
    ### make a return item when the Frag genus is not selected
    tax_string_i <- 'non-frag-genus'
  }
  return(tax_string_i)
}

get_tax_assoc_comprehensive <- function(quant_out, tfg_annot, tara_out){
  
  
  # quant_out <- all_avail_csvs_unique_only
  # tfg_annot <- annot_mcl
  # tara_out <- tara_tax
  ## function to assign taxonomy to peptides
  empty_tax_string <- c()
  # go through every peptide
  for(row_i in 1:nrow(quant_out)){
    
    # row_i <- 22171
    print(row_i)
    # if there is only one matched protein, take the taxonomic string from that protein
    if(quant_out$n_proteins[row_i] == 1){
      # first convert it to a character      
      protein_i <- as.character(quant_out[row_i, ]$protein)
      # if it was matched in the TARA set, then used the taxa string
      if(grepl(pattern = "TARA", x = protein_i)){
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        tax_string_i <- tara_out[tara_out$`Genome_Id final names` == protein_i_taxa, ]$Best_taxonomy_GENRE
        tax_i <- tax_string_i
      } else {
        # take the taxa stringfrom the TFG annotation file
        tax_i <- tfg_annot[tfg_annot$orf_id == protein_i, ]$grpnorm_taxgrp %>% as.character()
      }
      
      if(is.na(tax_i)){
        tax_i <- "no-tax-assignment"
      }
      
      empty_tax_string <- c(empty_tax_string, tax_i)
      
    }
    # if it matched more than one protein, the plot thickens.
    if(quant_out$n_proteins[row_i] > 1){
      # get the list of proteins
      protein_i <- parse_contigs(protein_quant_out = quant_out, row_i = row_i)
      
      ### first test if there are any TARA proteins in here, if not, summarize by grpnorm_taxgrp
      
      if(sum(grepl(pattern = "TARA", x = protein_i)) == 0){
        class_i <- classify_using_tfg(tfg_annot = tfg_annot, 
                                      protein_i = protein_i)
      }
      
      # this sum might look weird, it's because protein_i is a vector and the grepl statement takes a vector of TRUE/FALSE. The sum ofa 
      # vector of TRUE/FALSE is 1./0 respectively
      # if these proteins are all TARA matches
      if(sum(grepl(pattern = "TARA", x = protein_i)) == length(protein_i)){
        
        # taking the taxonomic string out of the protein_id
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        
        all_genres <- unique(tara_out[tara_out$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$Best_taxonomy_GENRE)
        if(length(all_genres) == 1){
          class_i <- all_genres
        } else if(length(all_genres) > 1){
          class_i <- 'mixture-of-taxonomic-annotations'
        }
      }
      # if there is a mixture of TARA and non-TARA sequences
      
      mixture_of_seqs <- (sum(grepl(pattern = "TARA", x = protein_i)) > 0) & 
        (sum(grepl(pattern = "TARA", x = protein_i)) < length(protein_i))
      if(mixture_of_seqs){
        # taking the taxonomic string out of the protein_id
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        
        all_genres <- unique(tara_out[tara_out$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$Best_taxonomy_GENRE)
        if(length(all_genres) == 0){
          ## somtimes the only ones that are hits are actually decoy sequences
          all_genres <- "decoy-seq"
        }
        ### we also want to look at the non-TARA sequences in this set of proteins
        protein_i_no_tara <- protein_i[!grepl("TARA", protein_i)]
        class_i_tfg_temp <- classify_using_tfg(tfg_annot = tfg_annot, 
                                      protein_i = protein_i_no_tara)
        ## if there are multiple genres from TARA
        if(length(all_genres) > 1){
          class_i <- 'mixture-of-taxonomic-annotations'
        } else if(class_i_tfg_temp == all_genres){
          class_i <- all_genres
        } else {
          class_i <- 'mixture-of-taxonomic-annotations'
        }
        
      }
        
      empty_tax_string <- c(empty_tax_string, class_i)
      if(length(empty_tax_string) != row_i){
        print('its breaking here')
      }
    }
  }
  return(empty_tax_string)
}

check_for_frag <- function(tax_i_string, ## Taxonomic string associated with protein
                           protein_string_i, ## Protein string
                           tfg_annotations = tfg_annot){
  ### this gets a frag species level designation
  if(tax_i_string == 'Fragilariopsis'){
    ### determine if it's from the TARA data
    if(grepl(pattern = "TARA", x = protein_string_i)){
      ## Get the specific genome name from the TARA dataset
      tax_out_sum <- sub("_[^_]+$", "", protein_string_i)
    } else {
    ### if not from TARA, then from TFG dataset
    ### this returns the whole taxonomic string associated
      tax_out_sum <- tfg_annotations[tfg_annotations$orf_id == protein_string_i, ]$best_tax_string %>% as.character()
    }
  } else {
    tax_out_sum <- 'non-frag-genus'
  }
  return(tax_out_sum)
}

get_tax_assoc_comprehensive_frag <- function(quant_out, tfg_annot, tara_out){
  
  # quant_out <- all_avail_csvs_unique_only
  # tfg_annot <- annot_mcl
  # tara_out <- tara_tax
  
  ## string to assign taxonomy to peptides
  empty_tax_string <- c()
  
  ## string to assign sub-genus level taxonomy to Frag-level genus peptides.
  empty_tax_string_frag <- c()
  
  quant_out_empty <- data.frame(peptide = character(),
                                protein = character(),
                                n_proteins = numeric(),
                                assigned_tax = character(),
                                assigned_frag = character())
  
  ## go through every peptide
  for(row_i in 1:nrow(quant_out)){
    
    ## append the peptide id
    temp_quant_df_peptide <- quant_out[row_i, ]$peptide
    ## this appends a list of proteins that match to the peptide
    temp_quant_df_protein <- quant_out[row_i, ]$protein
    ## this appends a number of proteins
    temp_quant_df_n_proteins <- quant_out[row_i, ]$n_proteins
    
    # row_i <- 22171
    print(row_i)
    
    # if there is only one matched protein, take the taxonomic string from that protein
    if(quant_out$n_proteins[row_i] == 1){
      
      # first convert it to a character      
      protein_i <- as.character(quant_out[row_i, ]$protein)
      
      # if it was matched in the TARA set, then used the taxa string
      if(grepl(pattern = "TARA", x = protein_i)){
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        tax_string_i <- tara_out[tara_out$`Genome_Id final names` == protein_i_taxa, ]$Best_taxonomy_GENRE
        tax_i <- tax_string_i
      } else {
        # take the taxa string from the TFG annotation file
        tax_i <- tfg_annot[tfg_annot$orf_id == protein_i, ]$grpnorm_taxgrp %>% as.character()
      }
      
      if(is.na(tax_i)){
        tax_i <- "no-tax-assignment"
        protein_i <- 'no-protein'
      }
      
      ## Get the Frag species, if it is indeed Frag genus.
      frag_group_i <- check_for_frag(tax_i_string = tax_i, 
                                   protein_string_i = protein_i, 
                                   tfg_annotations = tfg_annot)
      
      ## append the group level associations, and the Frag-species level associations.
      # empty_tax_string_frag <- c(empty_tax_string_frag, tax_i_frag)
      # empty_tax_string <- c(empty_tax_string, tax_i)
      
    }
    
    # if it matched more than one protein, the plot thickens.
    if(quant_out$n_proteins[row_i] > 1){
      # get the list of proteins
      protein_i <- parse_contigs(protein_quant_out = quant_out, row_i = row_i)
      
      ### first test if there are any TARA proteins in here, if not, summarize by grpnorm_taxgrp
      if(sum(grepl(pattern = "TARA", x = protein_i)) == 0){
        tax_i <- classify_using_tfg(tfg_annot = tfg_annot, 
                                      protein_i = protein_i)
        ### Get a summary of if all the proteins belong to one Frag species or several.
        frag_group_i <- classify_frag_using_tfg(tfg_annot = tfg_annot, 
                                                protein_i = protein_i, 
                                                class_i = tax_i)
      }
      
      # this sum might look weird, it's because protein_i is a vector and the grepl statement takes a vector of TRUE/FALSE. The sum of a 
      # vector of TRUE/FALSE is 1./0 respectively
      # if these proteins are all TARA matches
      if(sum(grepl(pattern = "TARA", x = protein_i)) == length(protein_i)){
        
        # taking the taxonomic string out of the protein_id
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        
        all_genres <- unique(tara_out[tara_out$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$Best_taxonomy_GENRE)
        if(length(all_genres) == 1){
          tax_i <- all_genres
        } else if(length(all_genres) > 1){
          tax_i <- 'mixture-of-taxonomic-annotations'
        }
        
        ## If it's a Frag-genus labelled protein, then get a string of all the genome IDs
        if(tax_i == 'Fragilariopsis'){
          frag_group_i <- paste(protein_i_taxa, 
                                collapse = '___')
        } else if(tax_i == 'mixture-of-taxonomic-annotations'){
          frag_group_i <- 'non-frag-genus'
        }
      }
      
      # if there is a mixture of TARA and non-TARA sequences
      mixture_of_seqs <- (sum(grepl(pattern = "TARA", x = protein_i)) > 0) & 
        (sum(grepl(pattern = "TARA", x = protein_i)) < length(protein_i))
      
      if(mixture_of_seqs){
        
        # taking the taxonomic string out of the protein_id
        protein_i_taxa <- sub("_[^_]+$", "", protein_i)
        
        all_genres <- unique(tara_out[tara_out$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$Best_taxonomy_GENRE)
        if(length(all_genres) == 0){
          ## somtimes the only ones that are hits are actually decoy sequences
          all_genres <- "decoy-seq"
        }
        
        ## get all TARA genomes of Frag
        
        ### we also want to look at the non-TARA sequences in this set of proteins
        protein_i_no_tara <- protein_i[!grepl("TARA", protein_i)]
        class_i_tfg_temp <- classify_using_tfg(tfg_annot = tfg_annot, 
                                               protein_i = protein_i_no_tara)
        
        ### get the Frag-species level classification for these non-TARA sequences
        frag_class_i_tfg_temp <- classify_frag_using_tfg(tfg_annot = tfg_annot, 
                                                         protein_i = protein_i_no_tara, 
                                                         class_i = class_i_tfg_temp)
        
        ## if there are multiple genres from TARA
        if(length(all_genres) > 1){
          tax_i <- 'mixture-of-taxonomic-annotations'
          frag_group_i <- 'non-frag-genus'
          
        } else if(class_i_tfg_temp == all_genres){
          tax_i <- all_genres
          
          ### get all the Frag genomes from the TARA set
          all_genomes <- unique(tara_out[tara_out$`Genome_Id final names` %in% as.character(protein_i_taxa), ]$`Genome_Id final names`)
          if(length(all_genomes) == 0){
            ## sometimes the only ones that are hits are actually decoy sequences
            all_genomes <- ""
            print('going in here')
          }
          
          frag_group_i <- paste(c(frag_class_i_tfg_temp, all_genomes), collapse = '___')
          
        } else {
          tax_i <- 'mixture-of-taxonomic-annotations'
          frag_group_i <- 'non-frag-genus'
        }
        
      }
      
      # empty_tax_string_frag <- c(empty_tax_string_frag, frag_group_i)
      # empty_tax_string <- c(empty_tax_string, class_i)
      # if(length(empty_tax_string) != row_i){
      #   print('its breaking here')
      # }
    }
    
    temp_quant_df <- data.frame(peptide = temp_quant_df_peptide,
                                protein = temp_quant_df_protein,
                                n_proteins = temp_quant_df_n_proteins,
                                assigned_tax = tax_i,
                                assigned_frag = frag_group_i)
    ## append df
    quant_out_empty <- rbind(quant_out_empty, temp_quant_df)
    
  }
  
  return(quant_out_empty)
}

###### ribosomal mass fraction

get_all_descriptions_for_contigs <- function(parsed_contigs_out, mcl_df){
  
  # collects all the descriptions associated with the contigs
  
  vector_of_desc_master <- c()
  
  vector_kegg <- c()
  vector_ko <- c()
  vector_kogg <- c()
  vector_pfams <- c()
  vector_tigrfams <- c()
  
  # now this function is greedy. as soon as it identifies a coarse grained match,
  # it assigns it to the protein coarse grain.
  for(j in 1:length(parsed_contigs_out)){
    # print(j)
    # j <- 1
    kegg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$kegg_desc
    ko_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KO_desc
    kogg_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$KOG_desc
    pfams_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$PFams_desc
    tigrframs_desc_contig <- mcl_df[mcl_df$orf_id == parsed_contigs_out[j], ]$TIGRFams_desc
    
    vector_of_desc <- c(kegg_desc_contig,
                        ko_desc_contig,
                        kogg_desc_contig,
                        pfams_desc_contig,
                        tigrframs_desc_contig)
    
    vector_kegg <- c(vector_kegg, kegg_desc_contig)
    vector_ko <- c(vector_ko, ko_desc_contig)
    vector_kogg <- c(vector_kogg, kogg_desc_contig)
    vector_pfams <- c(vector_pfams, pfams_desc_contig)
    vector_tigrfams <- c(vector_tigrfams, tigrframs_desc_contig)
    
    # print(temp_match)
    # print(unique(temp_match))
    # print(vector_of_desc)
    vector_of_desc_master <- c(vector_of_desc_master, vector_of_desc)
    # checking to see if there are any non-U values
  }
  
  desc_df_all <- data.frame(kegg_desc = paste(unique(vector_kegg), collapse = '____'),
                            ko_desc = paste(unique(vector_ko), collapse = '____'),
                            kogg_desc = paste(unique(vector_kogg), collapse = '____'),
                            pfams_desc = paste(unique(vector_pfams), collapse = '____'),
                            tigrfams_desc = paste(unique(vector_tigrfams), collapse = '____'))
  
  return(list(vector_of_desc_master, desc_df_all))
  
}

check_length_model_grains <- function(unique_coarse_grain){
  
  # intermediate function that takes the length of coarse grain assignments 
  # and processes them to determine which assignment it should be
  
  # if they are not all identical, assign the peptide to 'U'
  if(length(unique_coarse_grain) > 1){
    coarse_assignment <- "U"
  }
  # if they are all identical, then assign it
  if(length(unique_coarse_grain) == 1){
    coarse_assignment <- unique_coarse_grain
  }
  # if there are none, then assign it to a 0
  if(length(unique_coarse_grain) == 0){
    coarse_assignment <- "U"
  }
  
  return(coarse_assignment)
  
}

match_desc_to_model_grains <- function(description_vec_out){
  
  # matches description of output to model
  
  # if there are no descriptions, then return
  if(length(description_vec_out) == 0){
    coarse_assignment <- 'U'
  }
  
  model_coarse_grain <- c()
  
  if(length(description_vec_out) > 0){
    # matching descriptions to coarse grained models
    # goes through each description in the vector
    for(i in 1:length(description_vec_out)){
      
      ribosome_string <- c("ribosom")
      ribosome_anti_string <- c("synthesis")
      
      if(grepl(paste(ribosome_string, collapse = "|"), 
               description_vec_out[i], ignore.case = TRUE) &
         !grepl(paste(ribosome_anti_string, collapse = "|"), 
                description_vec_out[i], ignore.case = TRUE)){
        model_coarse_grain <- c(model_coarse_grain, "R")
      }
      
      # model_coarse_grain <- c("U", "R")
    }
    
    if(length(model_coarse_grain) == 0){
      model_coarse_grain <- 'U'
    }
    
    # check that all the coarse grains are identical
    unique_coarse_grain <- unique(model_coarse_grain)
    
    coarse_assignment <- check_length_model_grains(unique_coarse_grain = unique_coarse_grain)
    
  }
  
  return(coarse_assignment)
  
}

format_desc_df <- function(desc_out_intermediate){
  # need to add NA's if there are no annotations for any orfs returned  
  if(nrow(desc_out_intermediate) == 0){
    desc_out_intermediate[nrow(desc_out_intermediate)+1,] <- NA
  }
  return(desc_out_intermediate)
}

get_coarse_grains_ribo <- function(tax_assign_df, mcl_df){
  ### function to assign ribosomal proteins
  
  all_coarse_assignments <- c()
  
  all_associated_desc <- data.frame(kegg_desc = character(),
                                    ko_desc = character(),
                                    kogg_desc = character(),
                                    pfams_desc = character(),
                                    tigrfams_desc = character())
  
  for(pep_i in 1:nrow(tax_assign_df)){
    
    print(pep_i)
    # print(diatom_peps[i, ]$go_descs)
    # parse the output, getting all associated contigs associated with the peptide
    contigs_parsed <- parse_contigs(protein_quant_out = tax_assign_df,
                                    row_i = pep_i)
    # for each contig, get a vector of descriptions
    all_descs <- get_all_descriptions_for_contigs(parsed_contigs_out = contigs_parsed, 
                                                  mcl_df = mcl_df)
    
    # match descriptions to the model output
    coarse_assignment_val <- match_desc_to_model_grains(all_descs[[1]])
    
    desc_df <- format_desc_df(all_descs[[2]])
    
    # appending the description to the 
    all_associated_desc <- rbind(all_associated_desc, desc_df)
    all_coarse_assignments <- c(all_coarse_assignments, coarse_assignment_val)
    
  }
  
  tax_assign_df$ribosomal <- all_coarse_assignments
  tax_assign_df <- cbind(tax_assign_df, all_associated_desc)
  
  return(tax_assign_df)
}

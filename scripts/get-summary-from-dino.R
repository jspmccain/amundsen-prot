library(magrittr)

## first list all files from Dinosaur
# setwd('amundsen-prot/scripts/')

dino_files <- list.files(path = 'data/dino-converted/', pattern = "\\.tsv$")

## then loop through all the file names and get the total intensity sum for the injection
loop_through_files_get_int <- function(list_of_files){
  ## make aggregating dataframe
  out_df <- data.frame(file_name = character(),
                       sum_of_int = numeric())
  
  ## loop through
  for(i in 1:length(list_of_files)){
    ## add in the file path
    file_i <- paste0('data/dino-converted/', list_of_files[i])
    ## read in the file
    tsv_read <- read.table(file = file_i, header = TRUE)
    ## calculate sum intensity metric
    intensity_metric <- tsv_read$intensitySum %>% sum()
    
    ## aggregate this intermediate
    int_df <- data.frame(file_name = file_i,
                         sum_of_int = intensity_metric)
    ## append this intermediate
    out_df <- rbind(out_df, int_df)
  }
  
  return(out_df)
}

## use this function
dino_files_with_metric <- loop_through_files_get_int(list_of_files = dino_files)

write.csv(x = dino_files_with_metric, file = '../data/dino-converted/dino_summaries.csv', row.names = FALSE)

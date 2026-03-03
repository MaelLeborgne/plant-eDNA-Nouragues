
# ===== data_import.R =====
# Regroup importation functions by Subject

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # Basic
library(magrittr)  # %<>%
library(ROBIFastread) # Read fasta
library(Matrix) # output of OBITool
library(taxize) # taxid and NCBI
library(jsonlite) # manage json
library(vegan)

# ===== GENERAL =====
# Create reference dataframe for sample name and position
create_sample_ref <- function(PROJECT_PATH){
  # Import physical data of samples
  data_samples <- read.table(
    "/home/jean/Bureau/2026_Stage_ADNe/doi_10_5061_dryad_dc28m__v20181029/Zinger_et_al_metabar_PP1_samples.txt",
    header = TRUE,
    comment.char = "#"
  )
  # Select and split pcr
  sample_ref <- data_samples %>%
    select(pcr,x,y) %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    distinct(sample_loc, x, y)
  
  # Return Ref
  sample_ref
}

# ==== BOTA ====
import_bota <- function(PROJECT_PATH){
  # Import CSV
  bota <- read.csv(paste0(PROJECT_PATH,"data/Nour_2012.csv"),sep = ",",dec = ".") %>%
    mutate(Species = paste(Genus, Species, sep=" "))
  # Return
  bota
}

# ==== GENET ====
import_genet <- function(PROJECT_PATH){
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/"))
  
  # GH marker
  diet_data_GH <- read_obifasta("metabar_PP_GH.simple.clean.c10.l40.taxo.ann.fasta")
  diet_data_GH %<>% extract_features("obitag_bestmatch","obitag_bestid","obitag_rank","scientific_name",'taxid')
  
  # 18S marker
  diet_data_18S <- read_obifasta("metabar_PP_18s.simple.clean.c10.l80.taxo.ann.fasta")
  diet_data_18S %<>% extract_features("obitag_bestmatch","obitag_bestid","obitag_rank","scientific_name",'taxid')
  
  # ITS marker
  diet_data_ITS <- read_obifasta("metabar_PP_ITS.simple.clean.c10.l180.taxo.ann.fasta")
  diet_data_ITS %<>% extract_features("obitag_bestmatch","obitag_bestid","obitag_rank","scientific_name",'taxid')
  
  # Centralize data in Named List (no duplication, just ref)
  gene <- setNames(
    list(diet_data_GH, # df
         diet_data_18S, 
         diet_data_ITS), 
    c("GH","18S","ITS") # names
  )
  
  # Return
  gene
}

# Extract obiclean_weight (key)
extract_feature <- function(data, key){
  # with ROBIFastread function
  weight_mat <- extract_readcount(data, key=key)
  # The result is a dgCMatrix
  weight_mat
}

# Sum all different read of the focus taxa for each sample
sample_seq_count <- function(genet_focus){
  
  # Extract quantitative value in a matrix 
  sample_matrix <- extract_feature(genet_focus,"obiclean_weight")
  
  # Weight by sample
  sample_weight <- data.frame(
    sample_id = rownames(sample_matrix),
    Weight = rowSums(sample_matrix)
  )
  
  # Nombre d'OTU
  sample_otu <- data.frame(
    sample_id = rownames(sample_matrix),
    OTU = rowSums(sample_matrix > 0)
  )
  
  # Diversité de Shannon
  sample_shannon <- data.frame(
    sample_id = rownames(sample_matrix),
    Div_alpha = diversity(sample_matrix, index = "shannon")
  )
  
  # Join
  sample_df <- left_join(sample_otu, sample_shannon, by = "sample_id") %>%
    left_join(sample_weight, by = "sample_id")
  
  # Return the doublet (sample_id, Count)
  return(sample_df)
}

# Mean samples Count by localisation
loc_mean <- function(sample_df){
  
  # Mean samples count by localisation
  sample_df <- sample_df %>%
    separate(sample_id, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    group_by(sample_loc) %>%
    summarise(
      Weight    = mean(.data[["Weight"]]   , na.rm = TRUE),
      OTU       = mean(.data[["OTU"]]      , na.rm = TRUE),
      Div_alpha = mean(.data[["Div_alpha"]], na.rm = TRUE)
      )
  
  # Add (x,y) with sample_ref
  mean_sample <- left_join(sample_df, samples_ref, by="sample_loc")
  
  # Complete missing value with zeros
  genet_grid <- ref_grid() %>%
    left_join(mean_sample, by = c("x", "y")) %>%
    mutate(
      Weight       = replace_na(.data[["Weight"]],    0),
      OTU          = replace_na(.data[["OTU"]],       0),
      Div_alpha    = replace_na(.data[["Div_alpha"]], 0)
      )
  
  # Return
  return(genet_grid)
}


# ==== TAXID REF ====

# Function Take the rank name for a give rank
get_rank <- function(x, rank) {
  val <- x$name[x$rank == rank]
  if(length(val) == 0) NA else val
}

# Fonction to create taxid ref for local genet data from ncbi
create_taxid_ref <- function(){
  # Take all unique taxid from clean genet data
  all_unique_taxid <- c(gene$GH$taxid, 
                        gene$`18S`$taxid, 
                        gene$ITS$taxid) %>% unique()
  
  # from unique taxid import list of full path
  taxo <- classification(all_unique_taxid, "ncbi")
  
  # Import F, G, s from taxo thanks to get.rank()
  tax_table <- map_dfr(taxo, ~tibble(
    family  = get_rank(.x, "family"),
    genus   = get_rank(.x, "genus"),
    species = get_rank(.x, "species")
  ))
  
  # remove the genus from the specie
  tax_table$species <- word(tax_table$species, 2) # keep the second word
  
  # Add taxid to FGs
  tax_table$taxid <- all_unique_taxid
  
  # Return
  tax_table
}

# ==== CREATE AND SAVE DATA ====
# # Import functions
# bota <- import_bota(PROJECT_PATH) # Bota
# gene <- import_genet(PROJECT_PATH) # Genet
# tax_table <- create_taxid_ref() # Taxid
# samples_ref <- create_sample_ref(PROJECT_PATH) # Samples names / position
# 
# # Save in a file
# base::save.image("~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.DataImport")


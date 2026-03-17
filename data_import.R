
# ===== data_import.R =====
# Regroup importation functions by Subject

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # Basic
library(magrittr)  # %<>%
library(ROBIFastread) # Read fasta
library(Matrix) # output of OBITool
library(taxize) # taxid and NCBI
library(jsonlite) # manage json
library(vegan)

# ===== Source data prep =====
source("data_prep.R")

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
sample_seq_count <- function(genet_focus = select_genet()){
  
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
    Shannon = diversity(sample_matrix, index = "shannon")
  ) %>% mutate(H1 = exp(Shannon))
  
  # Join per samples
  sample_df <- left_join(sample_otu, sample_shannon, by = "sample_id") %>%
    left_join(sample_weight, by = "sample_id")
  
  
  # --- Résumé global ---
  otu_tab <- colSums(sample_matrix)
  total_div <- diversity(otu_tab,index="shannon")
  
  summary_df <- data.frame(
    sample_id = "TOTAL_",
    Shannon = round(total_div,2)
  )
  
  # Combiner la grille et le résumé
  grid_with_summary <- bind_rows(sample_df, summary_df)
  
  # Retourner la liste : per-sample + résumé global
  return(grid_with_summary)
}

# Mean samples Count by localisation
loc_mean <- function(sample_df = sample_seq_count()){
  
  # Mean samples count by localisation
  sample_df <- sample_df %>%
    separate(sample_id, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    group_by(sample_loc) %>%
    summarise(
      Weight    = mean(.data[["Weight"]]   , na.rm = TRUE),
      OTU       = mean(.data[["OTU"]]      , na.rm = TRUE),
      Shannon = mean(.data[["Shannon"]], na.rm = TRUE),
      H1 = mean(.data[["H1"]], na.rm = TRUE)
      )
  
  # Add (x,y) with sample_ref
  mean_sample <- left_join(sample_df, samples_ref, by="sample_loc")
  
  # Complete missing value with zeros
  genet_grid <- ref_grid() %>%
    left_join(mean_sample, by = c("sample_loc","x", "y")) %>%
    mutate(
      Weight       = replace_na(.data[["Weight"]],    0),
      OTU          = replace_na(.data[["OTU"]],       0),
      Shannon    = replace_na(.data[["Shannon"]], 0),
      H1    = replace_na(.data[["H1"]], 0)
      ) %>%
    mutate(
      Weight_density = normalization(Weight)
      ) %>% bind_rows(filter(sample_df, sample_loc =="TOTAL"))
  
  # Return
  return(genet_grid)
}


# ==== TAXID REF ====

# Function Take the rank name for a give rank
get_rank <- function(x, rank) {
  val <- x$name[x$rank == rank]
  if(length(val) == 0) NA else val
}

# Get classification
get_all_classification <- function(){
  # Take all unique taxid
  all_unique_taxid <- c(gene$GH$taxid,
                        gene$`18S`$taxid,
                        gene$ITS$taxid) %>% unique()
  
  # Get full taxonomy
  taxo <- classification(all_unique_taxid, "ncbi")
  return(taxo)
}

# Create taxid ref
create_taxid_ref <- function(taxo){
  
  # Loop on each taxid
  tax_table <- map_dfr(names(taxo), function(id){
    
    x <- taxo[[id]]
    
    tibble(
      taxid   = as.integer(id),
      id_rank = tail(x$rank, 1),
      id_accuracy = length(x$rank),
      family  = get_rank(x, "family"),
      genus   = get_rank(x, "genus"),
      species = get_rank(x, "species"),
      
      # Detect if Embryophyta is present in lineage
      is_plant = any(x$name == "Embryophyta"),
      
      # monocot clade
      is_monocot = any(x$name %in% c("Liliopsida", "monocotyledons"))
    )
  })
  
  # No family identification
  tax_table$not_na_family <- !is.na(tax_table$family)

  # # Monocot
  # not_trees <- c("Poaceae","Arecaceae","Araceae")

  # Species
  tax_table$not_na_sp <- !is.na(tax_table$species)

  # Valid species name
  tax_table$valid_sp_name <- tax_table$species %in% tax_valid_sp(tax_table)
  
  # Is in bota inventory
  tax_table$is_in_bota <- tax_table$species %in% bota$Species
  
  # # Growth form
  library(readxl)
  plant_growth <- read_excel("/home/jean/Bureau/ADNe_Nourague_2026/METABAR_plantes/metadata/Plants_Nouragues2001.xls")

  tax_table$exist_in_plant_growth <- tax_table$species %in% plant_growth$species

  return(tax_table)
}

test_tax_table_vec <- function(taxo, vec_names){
  
  tax_table <- map_dfr(names(taxo), function(id){
    
    x <- taxo[[id]]
    
    tibble(
      taxid   = as.integer(id),
      
      # monocot clade
      is_in_vec = any(x$name %in% vec_names)
    )
  })
  
  result <- table(tax_table$is_in_vec)
  
  return(result)
}

test_tax_table_2names <- function(taxo, name1 ="Arecaceae", name2="Liliopsida"){
  
  tax_table <- map_dfr(names(taxo), function(id){
    
    x <- taxo[[id]]
    
    tibble(
      taxid   = as.integer(id),
      
      # monocot clade
      is_in_vec = any(x$name == name1) & !any(x$name == name2)
    )
  })
  
  res_tab <- table(tax_table$is_in_vec)
  res_taxid <- tax_table %>% filter(is_in_vec) 
  
  return(list(res_tab,res_taxid))
}

check_gf <- function(sp){
  
  res <- occ_search(
    scientificName = sp,
    country = "GF",
    limit = 1
  )
  nrow(res$data) > 0
}

check_gf_all <- function(){
  library(rgbif)
  
  # Clean species name
  sp_clean <- tax_sp_vec(tax_table)
  
  # Get GBIF usage Keys
  gbif_match <- name_backbone(name = sp_clean)
  
  # Get GBIF count of observations for each species
  presence_GF <- map_int(
    gbif_match$usageKey,
    ~occ_count(taxonKey = .x, country = "GF")
  )
  # check if count > 0
  tax_table$in_guyane <- presence_GF > 0
  
  # Apply test
  gf_presence <- tibble(
    species = sp_clean,
    present_GF = sapply(sp_clean, check_gf)
  )
}

clean_species <- function(x){
  
  x %>%
    
    # enlever cf / aff
    gsub(" cf\\..*", "", .) %>%
    gsub(" aff\\..*", "", .) %>%
    
    # # enlever sp.
    # gsub(" sp\\..*", "", .) %>%
    
    # enlever collecteurs
    gsub(" [A-Z][a-z]+ [0-9].*", "", .) %>%
    
    # enlever hybrides
    gsub(" x .*", "", .) %>%
    
    trimws()
}

tax_valid_sp <- function(tax_table){
  sp_vec <- tax_table %>% filter(!is.na(species)) %>% pull(species)
  # Remove out format
  sp_clean <- clean_species(sp_vec)
  # Garder seulement les binomes
  sp_clean <- sp_clean[grepl("^[A-Z][a-z]+ (sp\\.|[a-z-]+)$", sp_clean)]
  return(sp_clean)
}


# # Plants zinger
# setwd(paste0(PROJECT_PATH,"zinger/"))
# zinger_plant <- read_obifasta("Zinger_et_al_metabar_PP1_Plants.uniq.fasta")

# ===== GENERAL =====

# Create reference dataframe for sample name and position
create_sample_ref <- function(PROJECT_PATH){
  # Import physical data of samples
  data_samples <- read.table(
    "/home/jean/Bureau/ADNe_Nourague_2026/zinger/Zinger_et_al_metabar_PP1_samples.txt",
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

# ==== ABIOTIC ====

import_abiotic <- function(){
  
  data_samples <- read.table(
    "/home/jean/Bureau/ADNe_Nourague_2026/zinger/Zinger_et_al_metabar_PP1_samples.txt",
    header = TRUE,
    comment.char = "#"
  )
  
  # Select and split pcr
  abiotic <- data_samples %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    filter(!control)
  
  # Return Ref
  return(abiotic)
}

# ==== CREATE AND SAVE DATA ====
# # Create
bota <- import_bota(PROJECT_PATH)   # Bota
gene <- import_genet(PROJECT_PATH)  # Genet
taxo <- get_all_classification()    # Classification genet
tax_table <- create_taxid_ref(taxo) # Taxid
samples_ref <- create_sample_ref()  # Samples names / position
abiotic <- import_abiotic()         # Abiotic

# List
save(
  bota,
  gene,
  taxo,
  tax_table,
  samples_ref,
  abiotic,
  file = "~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.DataImport.RData"
)


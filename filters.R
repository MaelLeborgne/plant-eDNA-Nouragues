
# ===== filters.R =====

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)     # Basic
library(readxl)    # Excel


# ===== TAXID REF =====

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
  
  # Species
  tax_table$not_na_sp <- !is.na(tax_table$species)
  
  # Valid species name
  tax_table$valid_sp_name <- tax_table$species %in% tax_valid_sp(tax_table)
  
  # Is in bota inventory
  tax_table$is_in_bota <- tax_table$species %in% bota$Species
  
  # Growth form
  plant_growth <- read_excel("/home/jean/Bureau/ADNe_Nourague_2026/METABAR_plantes/metadata/Plants_Nouragues2001.xls")
  
  tax_table$exist_in_plant_growth <- tax_table$species %in% plant_growth$species
  
  return(tax_table)
}

# ---- Test on Monocot remove ----

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

# ---- Species name format ----

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

# ---- Filtre Géographique ----

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



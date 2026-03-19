
# ===== data_prep.R =====
# Regroup data preparation functions by Subject

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)   # Basic

# ===== 1. FAMILY COMPOSITION BY MARKER =====

give_level_names <- function(config){
  
  if(config$Type == "genet"){
    # Import genet
    genet <- import_genet(config)
    
    # names in the level
    l <- tolower(config$Level) # tolower
    names <- unique(genet$motus[[l]]) # unique
    names <- names[!is.na(names)] # Remove NA
  }
  else if(config$Type == "bota"){
    # Import Bota
    bota <- import_bota(config)
    
    # names in the level
    names <- unique(bota[[config$Level]])
  }
  else{
    stop(paste0("Unknown config$Type : ", config$Type))
  }
  
  return(names)
}

filter_select_level <- function(Level){
  filters <- list()
  if(Level == "Family"){
    filters$filter_na_family = TRUE
  }
  else if(Level == "Genus"){
    filters$filter_na_family = TRUE
  }
  else if(Level == "Species"){
    filters$filter_na_family = TRUE
  }
  else{
    stop(paste0("Unknown Level: ", Level))
  }
  return(filters)
}

summary_ab_rel <- function(config){
  
  # Config All
  config_all <- config
  config_all$Name <- "All"
  filters <- filter_select_level(config$Level)
  
  grid_all <- grid_read(config_all, filters)
  total <- sum(grid_all$df$Metric)
  
  # Récupère tous les noms du clade
  all_names <- give_level_names(config)
  
  results <- lapply(all_names, function(n) {
    print(n)
    
    # Config Name
    config_name <- config
    config_name$Name <- n
    
    # Lit les données
    grid <- grid_read(config_name)
    
    # Total
    total_name <- sum(grid$df$Metric)
    ab_rel <- (total_name*100) / total
    
    data.frame(
      Clade = n,
      Metric = ab_rel
    )
  })
  
  final_df <- do.call(rbind, results)
  return(final_df)
}

# percent_other <- function(df, var_ref, var_per_cols, threshold){
#   
#   # lignes sous le seuil
#   df_low <- df %>%
#     dplyr::filter(.data[[var_per_cols[1]]] < threshold)
#   
#   # somme des colonnes de pourcentage
#   other_percent <- df_low %>%
#     dplyr::summarise(across(all_of(var_per_cols), sum, na.rm = TRUE))
#   
#   # ligne "Other"
#   line_other <- other_percent %>%
#     dplyr::mutate(!!var_ref := paste0("<", threshold, "%")) %>%
#     dplyr::select(all_of(var_ref), everything())
#   
#   # retirer lignes < threshold et ajouter Other
#   df <- df %>%
#     dplyr::filter(.data[[var_per_cols[1]]] >= threshold) %>%
#     dplyr::bind_rows(line_other)
#   
#   return(df)
# }
# 
# percent_out_bota <- function(df, 
#                              var_ref ="Clade", 
#                              var_per = "Genet_ab_rel", 
#                              var_na = "Bota_ab_rel"){
#   
#   # Somme des pourcentages hors bota
#   other_percent <- df %>%
#     dplyr::filter(is.na(.data[[var_na]])) %>%
#     dplyr::summarise(sum = sum(.data[[var_per]], na.rm = TRUE)) %>%
#     dplyr::pull(sum)
#   
#   # Ligne Out_of_bota
#   line_other <- tibble::tibble(
#     !!var_ref := "Out_of_bota",
#     !!var_per := other_percent
#   )
#   
#   # Conserver uniquement les lignes avec Bota
#   df <- df %>%
#     dplyr::filter(!is.na(.data[[var_na]])) %>%
#     dplyr::bind_rows(line_other)
#   
#   return(df)
# }

# ===== 2. CORRELATION TABLE =====

bota_names <- function(){
  
  # Import Bota
  config <- list(
    Year = 2012
  )
  bota <- import_bota(config)
  
  # List bota
  list_bota <- list(
    "Family" = unique(bota[["Family"]]),
    "Genus" = unique(bota[["Genus"]]),
    "Species" = unique(bota[["Species"]])
    )
  
  return(list_bota)
}

make_cor_table_bota_genet <- function(){
  
  # Make tab of configurations
  config_tab <- expand.grid(
    "g.Rar" = c(FALSE),
    "g.Marker" = c("GH","ITS","18S"),
    "g.Metric" = c("Reads"),
    "b.Metric" = c("Count"),#,"NCI"),
    "b.Radius" = c(10),
    "b.Year" = c(2012),
    "Level" = c("Family"),#, "Genus", "Species"),
    stringsAsFactors = FALSE
  )
  
  # Read "Name" for all Levels
  list_bota <- bota_names()
  # Add Name by Level
  config_tab <- config_tab %>%
    rowwise() %>%                        # travailler ligne par ligne
    mutate(Name = list(list_bota[[Level]])) %>%  # créer une liste de Names
    unnest(Name)                          # transformer la liste en lignes
  
  for(i in 1:nrow(config_tab)){
    print(i)
    # Compute grid correlation
    config_tab$Correlation[[i]] <- grid_correlation(config_tab[i,])
  }
  
  # length(list_bota$Family)
  # 63
  # Level = Family
  # 378
  # Level = All
  # 5346
  # TODO:
  # - Enlever les espèce/genus/family absente de génèt
  # - Enlever les famille de moin de 1% d'abondance relative (refaire la tab_ab_rel)
  
  # TEST
  # test_line <- config_tab[56,]
  # grid_correlation(test_line)
  # Erreur dans rowSums(reads_filtered) : 
  # 'x' doit être un tableau ayant au moins 2 dimensions
  
  return(config_tab)
}

grid_correlation <- function(row){
  
  # Extract configs
  configs <- extract_configs(row)
  
  # Read Bota grid
  grid_bota <- grid_read(configs$bota)
  
  # Read Genet grid
  grid_genet <- grid_read(configs$genet)
  
  # calcul grid correlation
  grid_corr <- cor(grid_bota$df$Metric, 
                   grid_genet$df$Metric)
  
  return(grid_corr)
  
}

extract_configs <- function(row) {
  
  # Colonnes botaniques : celles qui commencent par b. ou Level / Name
  b_cols <- grep("^b\\.|Level|Name", names(row), value = TRUE)
  bot_config <- as.list(row[, b_cols])
  bot_config$Type <- "bota"
  
  # Nettoyer les noms : enlever le préfixe b.
  names(bot_config) <- sub("^b\\.", "", names(bot_config))
  
  # Colonnes génétiques : celles qui commencent par g. ou Level / Name
  g_cols <- grep("^g\\.|Level|Name", names(row), value = TRUE)
  genet_config <- as.list(row[, g_cols])
  genet_config$Type <- "genet"
  
  # Nettoyer les noms : enlever le préfixe g.
  names(genet_config) <- sub("^g\\.", "", names(genet_config))
  
  list(
    bota  = bot_config,
    genet = genet_config
  )
}



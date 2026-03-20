
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
  # Read grid All
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

test_ab_rel_sum <- function(config_table){
  # Calcul tabs
  config_table$Result <- lapply(seq_len(nrow(config_table)), function(i) {
    message("Iteration ", i)
    summary_ab_rel(as.list(config_table[i, ]))
  })
  
  # Check 100%
  config_table$Check_100 <- lapply(seq_len(nrow(config_table)), function(i) {
    sum(config_table$Result[[i]]$Metric)
  })
  
  return(config_table)
}

percent_other <- function(df, var_ref, var_per_cols, threshold){

  # lignes sous le seuil
  df_low <- df %>%
    dplyr::filter(.data[[var_per_cols[1]]] < threshold)

  # somme des colonnes de pourcentage
  other_percent <- df_low %>%
    dplyr::summarise(across(all_of(var_per_cols), sum, na.rm = TRUE))

  # ligne "Other"
  line_other <- other_percent %>%
    dplyr::mutate(!!var_ref := paste0("<", threshold, "%")) %>%
    dplyr::select(all_of(var_ref), everything())

  # retirer lignes < threshold et ajouter Other
  df <- df %>%
    dplyr::filter(.data[[var_per_cols[1]]] >= threshold) %>%
    dplyr::bind_rows(line_other)

  return(df)
}

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


# Common Family for RAREFIED data
give_common_family <-function(){
  
  # Join result from ab_rel_tabs
  res1 <- inner_join(df_list[[1]],df_list[[7]],by="Clade")
  res3 <- inner_join(df_list[[3]],df_list[[7]],by="Clade")
  res5 <- inner_join(df_list[[5]],df_list[[7]],by="Clade")
  res_tot <- list(res1,res3,res5) %>% reduce(full_join, by = "Clade")
  
  # Take clade names
  names <- res_tot$Clade
  names <- names[names != "<2%"]
  
  return(names)
}

# UN-RAREFIED
give_genet_species <- function(Marker){
  
  # Config
  config <- list(
    Type = "genet",
    Rar = FALSE,
    Marker = Marker
  )
  # Import
  genet <- import_genet(config)
  
  # Genet not na species
  genet <- genet$motus %>%
    filter(!is.na(.data[["species"]])) # remove NA
  
  # Species vector
  spe_vec <- unique(genet[["species"]])
    
  return(spe_vec)
}

# Common species for UN-RAREFIED data
give_common_species <-function(){
  
  df_list <- list()
  
  # Import Bota
  config <- list(
    Year = 2012
  )
  df_list[[1]] <- import_bota(config) %>%
    filter(!is.na(Species)) %>%
    pull(Species)
  df_list[[1]] <- data.frame(Clade = unique(df_list[[1]]))
  
  # Import genet
  df_list[[2]] <- give_genet_species("GH")
  df_list[[2]] <- data.frame(Clade = df_list[[2]])
  df_list[[3]] <- give_genet_species("ITS")
  df_list[[3]] <- data.frame(Clade = df_list[[3]])
  df_list[[4]] <- give_genet_species("18S")
  df_list[[4]] <- data.frame(Clade = df_list[[4]])
  
  # Join result from ab_rel_tabs
  res2 <- inner_join(df_list[[1]],df_list[[2]],by="Clade")
  res3 <- inner_join(df_list[[1]],df_list[[3]],by="Clade")
  res4 <- inner_join(df_list[[1]],df_list[[4]],by="Clade")
  res_tot <- list(res2,res3,res4) %>% reduce(full_join, by = "Clade")
  
  # Take clade names
  names <- res_tot$Clade
  
  return(names)
}

make_cor_table_bota_genet <- function(config_tab, level_names){
  
  # Add Name by Level
  config_tab <- config_tab %>%
    rowwise() %>%                        # travailler ligne par ligne
    mutate(Name = list(level_names)) %>%  # créer une liste de Names
    unnest(Name)                          # transformer la liste en lignes
  
  # Init column
  config_tab$Correlation <- 0
  for(i in 1:nrow(config_tab)){
    print(i)
    # Compute grid correlation
    config_tab$Correlation[[i]] <- grid_correlation(config_tab[i,])
  }
  
  return(config_tab)
}

grid_correlation <- function(row){
  
  # Extract configs
  configs <- extract_configs(row)
  
  # Read Bota grid
  grid_bota <- grid_read(configs$bota)
  B <- grid_bota$df
  
  # Read Genet grid
  grid_genet <- grid_read(configs$genet)
  G <- grid_genet$df
  
  # Join
  BG <- inner_join(B, G, by = "sample_loc") 
  
  # True positive and true négative
  BG <- BG %>%
    mutate(tp = Metric.x == 1 & Metric.y == 1, # bota 1 genet 1
           tn = Metric.x == 0 & Metric.y == 0, # bota 0 genet 0
           )
  # Ratios
  Tp <- sum(BG$tp)/length(BG$tp)
  Tn <- sum(BG$tn)/length(BG$tn)
  
  if(!"Metric_spe" %in% names(row)){
    row$Metric_spe <- row$b.Metric
  }
  
  # Select Final metric
  grid_corr <- switch(row$Metric_spe,
                      "Count" = cor(B$Metric, G$Metric),
                      "NCI"   = cor(B$Metric, G$Metric),
                      "T+"    = Tp, # True positive
                      "T-"    = Tn, # True negative
                      "MCC"   = mcc(BG$Metric.x, BG$Metric.y),
                      stop(paste0("Unknown Metric :", configs$bota$Metric))
  )
  
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

mcc <- function(x, y) {
  TP <- sum(x == 1 & y == 1, na.rm = TRUE)
  TN <- sum(x == 0 & y == 0, na.rm = TRUE)
  FP <- sum(x == 0 & y == 1, na.rm = TRUE)
  FN <- sum(x == 1 & y == 0, na.rm = TRUE)
  
  # conversion en numérique pour éviter overflow
  TP <- as.numeric(TP)
  TN <- as.numeric(TN)
  FP <- as.numeric(FP)
  FN <- as.numeric(FN)
  
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  
  if (is.na(denom) || denom == 0) return(NA_real_)
  
  (TP * TN - FP * FN) / denom
}



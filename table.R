
# ===== data_prep.R =====
# Regroup data preparation functions by Subject

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)   # Basic

# ===== 1. TAXONOMY =====

# ---- Intersection table ----
extract_names <- function(df){
  
  cols <- c(
    "Family", "family",
    "Genus",  "genus",
    "Species","species"
  )
  
  col <- intersect(cols, names(df))[1]
  
  if(is.na(col)) {
    return(df)
  }
  
  df[[col]]
}

clean_names <- function(x){
  x %>%
    na.omit() %>%
    unique()
}

get_ab_tib <- function(level, marker){
  
  if(marker == "Bota"){
    
    df <- import_bota(list(Year = 2012))
    
    return(
      rank_ab_tib_bota(df, level)
    )
  }
  
  if(marker %in% c("GH", "ITS", "18S")){
    
    df <- import_genet(
      list(
        Marker = marker,
        Rar = TRUE
      )
    )
    
    return(
      rank_ab_tib_genet(df, tolower(level))
    )
  }
  
  if(marker == "REF"){
    
    return(
      import_plants_nouragues() %>%
        pull(tolower(level)) %>%
        unique()
    )
  }
  
  NULL
}

get_tax_list <- function(Level, Marker){
  # Bota
  if(Marker == "Bota"){
    df <- import_bota(list(Year = 2012))
    return(
      df %>% 
        pull(.data[[Level]]) %>%
        unique()
    )
  }
  # Genet
  if(Marker %in% c("GH", "ITS", "18S")){
    df <- import_genet(
      list(
        Marker = Marker,
        Rar = TRUE
      )
    )
    return(
      rank_ab_tib_genet(df, tolower(Level)) %>% 
        pull(.data[[tolower(Level)]])
    )
  }
  # Ref
  if(Marker == "REF"){
    return(
      import_plants_nouragues() %>%
        pull(tolower(Level)) %>%
        unique()
    )
  }
  
  # else
  NULL
  
}

# ---- Names and filters ----

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

# ---- Rank-Abundance Curve ----
# TIBLE
rank_ab_tib_bota <- function(df, Level){
  df %>% 
    group_by(.data[[Level]]) %>% 
    summarize(Abundance = n(), Dominance = n()/nrow(df)) %>%
    arrange(desc(Abundance)) %>%
    mutate(Rank = 1:length(Abundance))
}

rank_ab_tib_bota_BA <- function(df, Level){
  df %>% 
    group_by(.data[[Level]]) %>% 
    summarize(Abundance = sum(pi*(CircCorr/(2*pi))^2), Dominance = n()/nrow(df)) %>%
    arrange(desc(Abundance)) %>%
    mutate(Rank = 1:length(Abundance))
}

rank_ab_tib_genet <- function(df, Level){
  
  Reads <- colSums(df$reads)
  as.data.frame(Reads) -> Reads
  Reads$id <- row.names(Reads)
  
  left_join(df$motus, Reads, by = "id") %>%
    dplyr::filter(!is.na(.data[[Level]])) %>%
    group_by(.data[[Level]]) %>%
    summarise(Reads = sum(Reads))%>%
    ungroup() %>%
    mutate(Dominance = Reads / sum(Reads)) %>%
    arrange(desc(Reads)) %>%
    mutate(Rank = row_number())
}

# VECTOR
rank_ab_vec_bota <- function(Ab_tib, Level){
  
  # Build a named vector
  Ab_vec <- Ab_tib$Abundance
  names(Ab_vec) <- Ab_tib[[Level]]
  
  # Vector of abundances
  Ab_vec %>% as.AbdVector()
}

rank_ab_vec_genet <- function(Ab_tib, Level){
  # Build a named vector
  Ab_vec <- Ab_tib$Reads
  names(Ab_vec) <- Ab_tib[[Level]]
  
  # Vector of abundances
  Ab_vec %>% as.AbdVector()
}

# ---- Relative abondance ----

summary_ab_rel <- function(config){
  
  # Config All
  config_all <- config
  config_all$Name <- "All"
  
  # # Filter list
  # filters <- filter_select_level(config$Level) 
  # filters$filter_monocot <- FALSE
  
  # Read grid All
  grid_all <- grid_read(config_all)
  total <- sum(grid_all$df$Metric, na.rm = TRUE)
  
  # Récupère tous les noms du clade
  all_names <- c(give_level_names(config),"Indet.") %>% unique()
  
  results <- lapply(all_names, function(n) {
    print(n)
    
    # Config Name
    config_name <- config
    config_name$Name <- n
    
    # Lit les données
    grid <- grid_read(config_name)
    
    # Total
    total_name <- sum(grid$df$Metric, na.rm = TRUE)
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
    sum(config_table$Result[[i]]$Metric, na.rm = TRUE)
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

# ---- Growth form ----

growth_form_family <- function(fam_df){ 
  # for Plot
  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ .x / rowSums(across(starts_with("stat."))) * 100)) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat.")) %>%
    inner_join(fam_df, by="family")
}

tab_all_growth <- function(){
  
  # Return : all growth type table
  import_plants_nouragues() %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ round(.x / rowSums(across(starts_with("stat."))) * 100, digits = 2))) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat.")) %>%
    t() %>%
    as.data.frame() %>%
    arrange(desc(V1))
}

tab_family_growth <- function(){
  
  # Return : Table by family
  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ .x / rowSums(across(starts_with("stat."))) * 100)) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat.")) 
}


tab_family_ref <- function(){
  
  # Return : Nourague reference table for families %
  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE)))%>%
    mutate(
      total = rowSums(across(starts_with("stat.")), na.rm = TRUE)
    ) %>%
    dplyr::select(family,total) %>%
    ungroup() %>%
    mutate(
      total = (total / sum(total)) * 100
    ) %>%
    arrange(desc(total))
}

bias_assign <- function(df){
  df_long <- df %>%
    pivot_longer(
      cols = c(REF, GH, ITS, `18S`, BOTA),
      names_to = "inventory",
      values_to = "i"
    )
  
  ref <- df_long %>%
    dplyr::filter(inventory == "REF") %>%
    dplyr::select(Clade, i_ref = i)
  
  df_lr <- df_long %>%
    dplyr::filter(inventory != "REF") %>%
    dplyr::left_join(ref, by = "Clade") %>%
    dplyr::mutate(log_ratio = i/i_ref)#log((i + 1e-7) / (i_ref + 1e-7)))
  
  
  df_growth <- df_lr %>%
    pivot_longer(
      cols = c(terrestrial, tree, treelet, liana, shrub, vine, epilithic, parasite, hemi.epiphyte),
      names_to = "growth",
      values_to = "g"
    )
  
  df_bias <- df_growth %>%
    group_by(inventory, growth) %>%
    summarise(
      bias = weighted.mean(log_ratio, g, na.rm = TRUE)
    )
  
  return(df_bias)
}

# ===== 2. SPATIAL PATTERNS =====

# ---- Give names ----

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

give_genet_taxa <- function(Rar, Level, Marker){
  
  # Configuration
  genet_config <- list(
    Type   = "genet",
    Rar    = Rar,
    Marker = Marker,
    Level  = Level
  )
  
  # Importation
  genet <- import_genet(genet_config)$motus
  
  # Formatage
  level <- tolower(Level)
  
  # Filtrage 
  names <- genet %>%
    dplyr::filter(!is.na(.data[[level]])) %>%
    pull(.data[[level]]) %>%
    unique()
  
  return(names)
}

give_common_taxa <- function(Rar, Level) {
  
  # Bota vector direct
  bota <- import_bota(list(Year = 2012)) %>%
    dplyr::filter(!is.na(.data[[Level]])) %>%
    pull(.data[[Level]]) %>%
    unique()
  
  markers <- c("GH", "ITS", "18S")
  
  # Intersections avec chaque marker
  common_list <- lapply(markers, function(m) {
    intersect(bota, give_genet_taxa(Rar, Level, m))
  })
  
  # Union des résultats
  unique(unlist(common_list))
}

# ---- Correlation ----

make_cor_table_spatial <- function(config_tab, level_list){
    
  # Add Name by Level
  config_tab <- config_tab %>%
    rowwise() %>%
    mutate(Name = list(level_list[[Level]])) %>%  # Créer une liste de name
    unnest(Name)                        # transformer la liste en lignes
  
  # Initialisation to zero
  config_tab$Correlation <- 0 
  
  # For each row in config tab
  for(i in 1:nrow(config_tab)){
    print(i)
    # Compute grid correlation
    config_tab$Correlation[[i]] <- grid_correlation(config_tab[i,])
  }
  
  return(config_tab)
}

grid_correlation <- function(row){
  
  # Split Metric_corr cases
  if(row$Metric_corr %in% c("Count","NCI")){
    row$b.Metric <- row$Metric_corr  # Bota
    row$g.Metric <- "Reads"          # Génèt
  }
  else if(row$Metric_corr %in% c("T+","T-","F+","F-","MCC","Sen","Spe")){
    row$b.Metric <- "Presence"   # Bota
    row$g.Metric <- "Presence"   # Génèt
  }
  
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
           fp = Metric.x == 0 & Metric.y == 1, # bota 0 genet 1
           fn = Metric.x == 0 & Metric.y == 1  # bota 1 genet 0
           )
  
  # Ratios
  Tp <- sum(BG$tp)/length(BG$tp)
  Tn <- sum(BG$tn)/length(BG$tn)
  Fp <- sum(BG$fp)/length(BG$fp)
  Fn <- sum(BG$fn)/length(BG$fn)
  
  # Select Final metric
  grid_corr <- switch(row$Metric_corr,
                      "Count" = cor(B$Metric, G$Metric, method="spearman"),
                      "NCI"   = cor(B$Metric, G$Metric, method="spearman"),
                      "OTU"   = cor(B$Metric, G$Metric, method="spearman"),
                      "T+"    = Tp, # True positive
                      "T-"    = Tn, # True negative
                      "F+"    = Fp, # True positive
                      "F-"    = Fn, # True negative
                      "Sen"   = Tp / (Tp + Fn), # Sensitivity
                      "Spe"   = Tn / (Tn + Fp), # Specificity
                      "MCC"   = mcc(BG$Metric.x, BG$Metric.y),
                      stop(paste0("Unknown Metric :", configs$bota$Metric))
  )
  
  return(grid_corr)
}

# ---- Other ----
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

# ==== 3. DIVERSITY PATTERNS ====

make_cor_table_diversity <- function(config_tab){
  
  # Initialisation to zero
  config_tab$Correlation <- 0 
  
  # For each row in config tab
  for(i in 1:nrow(config_tab)){
    print(i)
    # Compute grid correlation
    config_tab$Correlation[[i]] <- grid_correlation(config_tab[i,])
  }
  
  return(config_tab)
}


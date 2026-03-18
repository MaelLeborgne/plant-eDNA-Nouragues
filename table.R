
# ===== data_prep.R =====
# Regroup data preparation functions by Subject

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)   # Basic


# ==== Correlation ====

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



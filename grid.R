
# ===== grid.R =====

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)   # Basic
library(tibble)  # Tibble class
library(stringr) # Work on string
library(purrr)   # pmap_int()
library(vegan)   # OTU diversity
library(terra)

# ===== INIT GRID =====

# Initialize empty grid
ref_grid <- function(){
  
  # Create all positions
  df <- expand.grid(
    x = seq(10, 290, 10), # X position of samples
    y = seq(10, 390, 10)  # Y position of samples
  )
  
  # Import samples name reference
  samples_ref <- create_sample_ref()
  
  # Join samples name
  df <- left_join(df, samples_ref, by = c("x","y"))
  
  # Init Metric to zero
  df$Metric <- 0
  
  return(df)
}

init_grid <- function(config, filters = list()){
  
  # Empty list
  grid <- list()
  
  # Attach empty grid
  grid$df <- ref_grid()
  
  # Attach config
  grid$config <- config
  
  # Attach filters
  grid$filters <- filters
  
  return(grid)
}

grid_normalize <- function(df){
  df$Metric <- df$Metric / sum(df$Metric)
  return(df)
}

grid_logaritmize <- function(df){
  df$Metric <- log(df$Metric+1)
  return(df)
}

grid_filtering <- function(df){
  
  library(terra)
  
  # 1. xtabs
  mat <- xtabs(Metric ~ x + y, data = df)
  
  # 2. CONVERSION FORCÉE (clé)
  mat <- as.matrix(mat)
  mat <- unclass(mat)
  
  # 3. raster
  mat <- t(mat)
  r <- terra::rast(mat)
  terra::ext(r) <- c(
    min(df$x), max(df$x),
    min(df$y), max(df$y)
  )
  
  # 4. filtre médian
  r_filt <- terra::focal(
    r,
    w = matrix(1, 3, 3),
    fun = median,
    na.rm = TRUE
  )
  
  # 5. retour data.frame
  df_filt <- as.data.frame(r_filt, xy = TRUE)
  names(df_filt)[3] <- "Metric"
  
  # 6. remettre sample_loc
  df_filt <- merge(
    df_filt,
    df[, c("x", "y", "sample_loc")],
    by = c("x", "y"),
    all.x = TRUE
  )
  
  return(df_filt)
}

# ==== TAXA FILTERS ====
select_tax <- function(Data, Level, Name){
  
  if(Name == "All"){
    focus <- Data
    return(focus)
  }
  else{
    # Select the focus in Data
    focus <- Data[!is.na(Data[[Level]]) & Data[[Level]] == Name,]
    return(focus)
  }

}

apply_filters <- function(df,
                          # filter_plant = TRUE,
                          filter_monocot = FALSE,
                          filter_na_family = FALSE,
                          filter_na_genus = FALSE,
                          filter_na_sp = FALSE
                          # filter_valid_sp = TRUE
                          ){
  
  # if(filter_plant){
  #   df <- df %>% dplyr::filter(is_plant)
  # }
  # 
  # if(filter_monocot){
  #   df <- df %>% dplyr::filter(!is_monocot)
  # }
  
  if(filter_na_family){
    df <- df %>% dplyr::filter(not_na_family)
  }
  
  if(filter_na_genus){
    df <- df %>% dplyr::filter(not_na_genus)
  }
  
  if(filter_na_sp){
    df <- df %>% dplyr::filter(not_na_sp)
  }
  
  # if(filter_valid_sp){
  #   df <- df %>% dplyr::filter(valid_sp_name)
  # }
  
  return(df)
}

# ===== BOTA =====

bota_build_matrix <- function(G, B) {
  
  df <- G$df
  C  <- G$config
  R  <- C$Radius
  
  all_levels <- unique(B[[C$Level]])
  
  # matrice vide
  mat <- matrix(0,
                nrow = nrow(df),
                ncol = length(all_levels))
  
  colnames(mat) <- all_levels
  
  for (pixel in 1:nrow(df)) {
    
    px <- df$x[pixel]
    py <- df$y[pixel]
    
    dij <- sqrt((px - B$Xfield)^2 + 
                  (py - B$Yfield)^2)
    
    if (R == "pixel") {
      sel <- which(B$Xfield > px-5 &
                     B$Xfield < px+5 &
                     B$Yfield > py-5 &
                     B$Yfield < py+5)
    } else {
      sel <- which(dij < R & dij > 0)
    }
    
    if (length(sel) > 0) {
      counts <- table(B[[C$Level]][sel])
      mat[pixel, names(counts)] <- counts
    }
  }
  row.names(mat) <- df$sample_loc
  
  return(mat)
}

get_neighbors <- function(px, py, B, R){
  
  dij <- sqrt((px - B$Xfield)^2 +
                (py - B$Yfield)^2)
  
  dij[dij > 100] <- 100
  dij[dij < 0.1] <- 0.1
  
  if (R == "pixel") {
    sel <- which(B$Xfield > px-5 &
                   B$Xfield < px+5 &
                   B$Yfield > py-5 &
                   B$Yfield < py+5)
  } else {
    sel <- which(dij < R & dij > 0)
  }
  
  list(sel = sel, dij = dij)
}

bota_grid_read <- function(grid, filters = list()) {
  
  # Configuration
  G <- grid
  C <- G$config
  M <- C$Metric
  R <- C$Radius
  
  # Importation
  bota <- import_bota(C)
  
  # Replace NA by Indet. (same as bota survey)
  for (l in c("Family","Genus","Species")) {
    bota[[l]][is.na(bota[[l]])] <- "Indet."
  }
  
  # Selection
  B <- select_tax(bota, C$Level, C$Name)
  
  # # Filter
  # B <- rlang::exec(apply_filters, B, !!!filters)
  
  # Build Matrix
  mat <- bota_build_matrix(G, B)
  
  # Fill Metric
  grid$df$Metric <- switch(M,
                           "Count"    = rowSums(mat),
                           "Presence" = as.integer(rowSums(mat) > 0),
                           "OTU"      = rowSums(mat > 0),
                           "Shannon"  = vegan::diversity(mat, index="shannon"),
                           "H1"       = exp(vegan::diversity(mat, index="shannon")),
                           "HYP"      = {
                             df <- G$df
                             
                             sapply(1:nrow(df), function(pixel){
                               
                               px <- df$x[pixel]
                               py <- df$y[pixel]
                               
                               res <- get_neighbors(px, py, B, R=10000)
                               sel <- res$sel
                               dij <- res$dij
                               
                               if (length(sel) == 0) return(0)
                               
                               sum(((B$CircCorr[sel]/pi)^C$Alpha) / dij[sel],
                                   na.rm = TRUE)
                             })
                           },
                           "EXP" = {
                             df <- G$df
                             
                             sapply(1:nrow(df), function(pixel){
                               
                               px <- df$x[pixel]
                               py <- df$y[pixel]
                               
                               res <- get_neighbors(px, py, B, R)
                               sel <- res$sel
                               dij <- res$dij
                               
                               if (length(sel) == 0) return(0)
                               
                               sum(((B$CircCorr[sel]/pi)^C$Alpha) *
                                     exp(-C$Beta * dij[sel]),
                                   na.rm = TRUE)
                             })
                           },
                           stop("Unknown Metric: ", M)
  )
  
  # Logaritmization
  if(isTRUE(C$Log)){
    grid$df <- grid$df %>% grid_logaritmize()
  }
  
  # Filtrage
  if(isTRUE(C$Filt)){
    grid$df <- grid$df %>% grid_filtering()
  }
  
  # Normalization
  if(isTRUE(C$Norm)){
    grid$df <- grid$df %>% grid_normalize()
  }
  
  # Return
  return(grid)
}

# ==== GENET READ ====

genet_grid_read <- function(grid, filters = list()){
  
  # Configuration
  G <- grid
  C <- G$config
  M <- C$Metric
  
  # Import genet
  genet <- import_genet(C)
  
  # Replace NA by Indet. (same as bota survey)
  for (l in c("family","genus","species")) {
    genet$motus[[l]][is.na(genet$motus[[l]])] <- "Indet."
  }
  
  # Select OTU
  otu_list <- select_tax(genet$motus, tolower(C$Level), C$Name)
  
  # # Apply otu filters
  # otu_list <- rlang::exec(apply_filters, otu_list, !!!filters) 
  
  # get otu_id
  otu_list <- otu_list %>% filter(!is.na(id)) %>% pull(id)
  
  # Replace : and - by . # (18S anomaly)
  otu_list <- gsub(":", ".", otu_list)
  otu_list <- gsub("-", ".", otu_list)
  
  # Select Reads from OTU_id
  reads_filtered <- genet$reads[, colnames(genet$reads) %in% otu_list, drop = FALSE]
  
  # Init
  reads_tax <- reads_filtered
  # SEULEMENT POUR SHANNON ET H1 (not species levels)
  # Grouper par OTU de même C$Level == C$Name
  if(M %in% c("Shannon","H1") & C$Level != "Species"){
    
    # Créer le vecteur d'espèces aligné avec les lignes
    tax_names <- genet$motus[[tolower(C$Level)]][
      match(colnames(reads_filtered), genet$motus$id)
    ]
    # Gestion des NA
    valid <- !is.na(tax_names)
    
    # Agréger par Name
    reads_tax <- t(rowsum(
      t(reads_filtered[, valid]),
      group = tax_names[valid]
    ))
  }
  
  # Create dataframe for calcul
  metric_df <- data.frame(sample_loc = rownames(reads_tax))
  
  # Compute the metric
  metric_df$Metric <- switch(M,
                    "Reads"         = rowSums(reads_tax),
                    "OTU"           = rowSums(reads_tax > 0),
                    "Presence"      = as.integer(rowSums(reads_tax > 0) > 0),
                    "Shannon"       = diversity(reads_tax, index = "shannon"),
                    "H1"            = exp(diversity(reads_tax, index = "shannon")),
                    stop("Unknown Metric: ", M)  # gestion d'erreur
  )
  
  # Join by sample_loc
  grid$df <- left_join(grid$df, metric_df, by="sample_loc")
  
  # Merge Metric cols
  grid$df <- grid$df %>%
    mutate(Metric = coalesce(Metric.y, Metric.x)) %>%
    dplyr::select(-Metric.x, -Metric.y)
  
  # Logaritmization
  if(isTRUE(C$Log)){
    grid$df <- grid$df %>% grid_logaritmize()
  }
  
  # Filtrage
  if(isTRUE(C$Filt)){
    grid$df <- grid$df %>% grid_filtering()
  }
  
  # Normalization
  if(isTRUE(C$Norm)){
    grid$df <- grid$df %>% grid_normalize()
  }
  
  return(grid)
}

# ==== ABIOTIC ====

abiotic_grid_read <- function(grid){
  
  grid$df <- import_abiotic() %>% 
    dplyr::select(x,y,.data[[grid$config$Metric]]) %>%
    rename(Metric = grid$config$Metric)
  
  return(grid)
}


# ==== GRID READ ====

grid_read <- function(config, filters=list()){
  
  # Initialize Grid
  grid <- init_grid(config)
  
  # Fill Grid
  if(config$Type == "bota"){
    grid <- bota_grid_read(grid, filters)
  }
  else if (config$Type == "genet"){
    grid <- genet_grid_read(grid, filters)
  }
  else if (config$Type == "abiotic"){
    grid <- abiotic_grid_read(grid)
  }
  else{
    stop(paste0("Unknown config$Type : ", config$Type))
  }
  
  return(grid)
}

build_grids <- function(config_tab){
  
  # Iterate on rows
  grid_list <- lapply(seq_len(nrow(config_tab)), function(i){
    # get parameters
    config <- as.list(config_tab[i, ])
    # Import grid
    grid_read(config)
  })
  
  return(grid_list)
}

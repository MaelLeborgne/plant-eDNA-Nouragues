
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

# ==== TAXA FILTERS ====
select_tax <- function(Data, Level, Name){
  
  focus <- NULL
  
  if(Name == "All"){
    focus <- Data
  }
  else{
    # Select the focus in Data
    focus <- Data[Data[Level] == Name,]
  }
  
  # Return
  return(focus)
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
  
  return(mat)
}

bota_grid_read <- function(grid, filters = list()) {
  
  G <- grid
  C <- G$config
  M <- C$Metric
  
  bota <- import_bota(C)
  B <- select_tax(bota, C$Level, C$Name)
  B <- rlang::exec(apply_filters, B, !!!filters)
  
  mat <- bota_build_matrix(G, B)
  
  grid$df$Metric <- switch(M,
                           "Count"    = rowSums(mat),
                           "Presence" = as.integer(rowSums(mat) > 0),
                           "OTU"      = rowSums(mat > 0),
                           "Shannon"  = vegan::diversity(mat, index="shannon"),
                           "H1"       = exp(vegan::diversity(mat, index="shannon")),
                           stop("Unknown Metric: ", M)
  )
  
  return(grid)
}

# # Gaussian
# bota_grid_Gauss <- function(bota_focus, BETA, GAMMA, R = 20, user_norm=FALSE){
#   # Give calcul to generic function
#   bota_grid_engine(
#     bota_focus,
#     R = R,
#     var_name = "Gaussian",
#     normalize = user_norm,
#     fun_pixel = function(df, sel, dij){
#       sum(
#         (pi*(df$CircCorr[sel]/2)^2) *
#           exp(-BETA * dij[sel]^GAMMA),
#         na.rm = TRUE
#       )
#     }
#   )
# }
# 
# ### Fonction to minimize (gauss)
# rss_function <- function(par, ADNe_obs) {
#   # Prameters
#   beta <- par[1]
#   gamma <- par[2]
# 
#   # Prediction
#   pred_grid <- bota_grid_Gauss(beta, gamma)
# 
#   # Residuals
#   resid <- ADNe_obs - pred_grid$density
# 
#   # RSS
#   print(sum(resid^2))
#   print(beta)
#   sum(resid^2)
#   }
# 
# ### Launch Optimization (gauss)
# optim_gauss <- function(genet_grid) {
#   res <- optim(
#     par = c(0.2,1),
#     fn = rss_function,
#     ADNe_obs = genet_grid$density,
#     method = "L-BFGS-B",
#     lower = c(0,0),  # logical constraint (biological)
#   )
#   return(res)
# }

# bota_grid_OTU_matrice <- function(bota_focus, R = 15){
#   # Give calcul to generic function
#   bota_grid_engine(
#     bota_focus,
#     R = R,
#     var_name = "OTU",
#     normalize = FALSE,
#     return_comm = TRUE, # Matrice
#     fun_pixel = function(df, sel, dij){
#       table(df$Species[sel])
#     }
#   )
# }




# ==== GENET READ ====

genet_grid_read <- function(grid, filters = list()){
  
  # Configuration
  G <- grid
  C <- G$config
  M <- C$Metric
  
  # Import genet
  genet <- import_genet(C)
  
  # Select OTU
  otu_list <- select_tax(genet$motus, tolower(C$Level), C$Name)
  # Apply otu filters
  otu_list <- rlang::exec(apply_filters, otu_list, !!!filters) 
  # get otu_id
  otu_list <- otu_list %>% filter(!is.na(id)) %>% pull(id)     
  # Replace : and - by . # (18S anomaly)
  otu_list <- gsub(":", ".", otu_list)
  otu_list <- gsub("-", ".", otu_list)
  
  # Select Reads from OTU_id
  reads_filtered <- genet$reads[, colnames(genet$reads) %in% otu_list, drop = FALSE]
  
  # Grouper par OTU de même C$Level == C$Name
  # Créer le vecteur d'espèces aligné avec les lignes
  tax_names <- genet$motus[[tolower(C$Level)]][
    match(colnames(reads_filtered), genet$motus$id)
  ]
  # Gestion des NA
  valid <- !is.na(tax_names)
  # Agréger par espèce
  reads_tax <- t(rowsum(
    t(reads_filtered[, valid]),
    group = tax_names[valid]
  ))
  
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


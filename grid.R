
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

init_grid <- function(config, filters){
  
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

# ==== SELECT TAXA ====
select_tax <- function(Data, Level, Name){
  
  focus <- NULL
  
  if(Name=="All"){
    focus <- Data
  }
  else{
    # Select the focus in Data
    focus <- Data[Data[Level] == Name,]
  }
  
  # Return
  return(focus)
}

# ===== BOTA =====

# Generic function to compute bota grid
bota_grid_engine <- function(B, G, fun_pixel
                             ){
  
  # Configuration
  df <- G$df
  C  <- G$config
  R  <- C$Radius
  
  # Pixel Loop
  for (pixel in 1:nrow(df)) {
    
    # adresse du pixel
    px <- df$x[pixel]
    py <- df$y[pixel]
    
    # distance de tout les arbre avec le pixel
    dij <- sqrt((px - B$Xfield)^2 + 
                  (py - B$Yfield)^2)
    
    # Sélectionner les arbres dans un intervale de distance
    sel <- which(dij < R & dij > 0)
    
    # fill the pixel with fun_pixel
    if (length(sel) != 0) {
        df$Metric[pixel] <- fun_pixel(B, sel, dij)
    }
  }
  return(df)
}

# Apply conditions for bota grid calcul
bota_grid_read <- function(grid){
  
  # Configuration
  G <- grid
  C <- G$config
  R <- C$Radius
  M <- C$Metric
  
  # Import bota
  bota <- import_bota(C)
  
  # Select Taxa
  B <- select_tax(bota, C$Level, C$Name)
  
  # Choix de la métrique botanique via switch
  grid$df <- switch(M,
               "Count"         = bota_grid_engine(B, G,
                                                  function(df, sel, dij){
                                                    length(sel)
                                                    }),
               "NCI"           = bota_grid_engine(B, G,
                                                  function(df, sel, dij){
                                                    sum((df$CircCorr[sel]^2) * 
                                                          exp(-0.25 * dij[sel]),
                                                        na.rm = TRUE)
                                                    }),
               "Basal_Area_m2" = bota_grid_engine(B, G,
                                                  function(df, sel, dij){
                                                    sum((bota[sel, ]$CircCorr^2) / (4 * pi * 10000))
                                                    }),
               "OTU"           = bota_grid_engine(B, G, 
                                                  function(df, sel, dij){
                                                    length(unique(df$Species[sel]))
                                                    }),
               "Shannon"       = bota_grid_engine(B, G,
                                               function(df, sel, dij){
                                                 sp_counts <- table(df$Species[sel])
                                                 vegan::diversity(sp_counts, index = "shannon")
                                                 }),
               "H1"            = bota_grid_engine(B, G,
                                              function(df, sel, dij){
                                                sp_counts <- table(df$Species[sel])
                                                exp(vegan::diversity(sp_counts, index = "shannon"))
                                              }),
               # "Gaussian"  = {
               #   result <- optim_gauss(B)
               #   bota_grid_Gauss(B, BETA = result$par[1], GAMMA = result$par[2])
               # },
               stop("Unknown var_bota: ", M)  # gestion d'erreur
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


# ===== GENET =====

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

# Get Familly, Genus, species from faxid reference
get_classification <-function(df){
  df <- inner_join(df, tax_table, by="taxid") 
  return(df)
}

# Apply filters
apply_filters <- function(df,
                          filter_plant = TRUE,
                          filter_monocot = TRUE,
                          filter_na_family = TRUE,
                          filter_na_sp = TRUE,
                          filter_valid_sp = TRUE){
  
  if(filter_plant){
    df <- df %>% dplyr::filter(is_plant)
  }
  
  if(filter_monocot){
    df <- df %>% dplyr::filter(!is_monocot)
  }
  
  if(filter_na_family){
    df <- df %>% dplyr::filter(not_na_family)
  }
  
  if(filter_na_sp){
    df <- df %>% dplyr::filter(not_na_sp)
  }
  
  if(filter_valid_sp){
    df <- df %>% dplyr::filter(valid_sp_name)
  }
  
  return(df)
}

# Select genetical reads with Marker and Name at a give Level
select_genet <- function(genet_all = gene, 
                         Marker = "GH", 
                         level = "family", 
                         Name = "All",
                         filters = list()){
  
  genet_focus <- genet_all[[Marker]] %>%
    get_classification()
  
  # appliquer filtres seulement si présents
  genet_focus <- rlang::exec(apply_filters, genet_focus, !!!filters)
  
  genet_focus <- genet_focus %>%
    select_tax(level, Name) %>%
    dplyr::filter(!is.na(.data[[level]]))
  
  return(genet_focus)
}

# Create genet grid for quantitative try
genet_grid_count <- function(genet_focus = select_genet()){
  
  # Create Grid
  grid <- ref_grid()
  
  # Get count of sequences by samples
  sample_df <- sample_seq_count(genet_focus)
  
  # Mean sample count by localisation
  genet_grid <- loc_mean(sample_df)
  
  # Remove NA sample_loc
  genet_grid <- genet_grid %>%
    filter(!is.na(sample_loc))
  
  # Add informations to TOTAL line
  genet_grid <- genet_grid %>% 
    mutate(
      Weight = ifelse(
        sample_loc == "TOTAL",
        sum(Weight[sample_loc != "TOTAL"]),
        Weight
      )
    )
  
  # Return
  return(genet_grid)
}

genet_grid_write <- function(focus, fileName){
  # Make grid
  grid <- genet_grid_count(focus)
  
  # Write
  write.csv(grid, paste0(PROJECT_PATH,"genet_grid/", fileName))
}

give_level_names <- function(Marker, Level, filters){
  level <- tolower(Level)
  # Focus marker
  mark_focus <- gene[[Marker]] %>% # Select the marker
    get_classification() %>% # Get Familly, Gender, species cols
    filter(!is.na(.data[[level]])) # remove NA
  
  mark_focus <- rlang::exec(apply_filters, mark_focus, !!!filters)
  # names in the level
  return(unique(mark_focus[[level]]))
}

genet_grid_write_all <- function(genet_all = gene, 
                                 markers  ="GH", 
                                 m_levels = "Family"){
  # For each markers
  for (m in markers) {
    # Focus marker
    mark_focus <- genet_all[[m]] %>% # Select the marker
      get_classification() # Get Familly, Gender, species cols
    # For each levels
    m_levels <- tolower(m_levels)
    for (l in m_levels) {
      L <- str_to_sentence(l) # upper first letter
      # remove NA
      mark_focus <- mark_focus %>%
        filter(!is.na(.data[[l]])) 
      # names in the level
      l_names <- c("All",unique(mark_focus[[l]]))
      # For each names
      for (n in l_names) {
        N <- str_to_sentence(n) # upper first letter
        # File name
        fileName <- paste0(m,"_",L,"_",N,".csv")
        print(fileName)
        # Focus taxa
        focus <- data.frame()
        if (N == "All"){
          focus <- mark_focus
        }
        else{
          focus <- select_tax(mark_focus, l, n) %>%
            filter(!is.na(.data[[l]])) # col family, genus or specie
        }
        # Write
        genet_grid_write(focus, fileName)
      }
    }
  }
}

# ==== WRITE ALL GENET_GRID ====
# # Test loop
# markers <- c("GH","18S","ITS")
# m_levels <- c("Family","Genus","Species")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)

# Test write
# markers <- c("GH")
# m_levels <- c("Family")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)

# ==== GENET READ ====

read_genet_grid_csv <- function(path, m, l, n){
  fileName <- paste0(m,"_",l,"_",n,".csv")
  read.csv(file.path(path, "genet_grid_5filters", fileName))
}

load_edna_data <- function(rarefied){
  if(!rarefied) return(NULL)

  dna_comm <- readRDS("...rds")
  dna_comm$reads
}

build_genet_grid <- function(m, l, n, filters){
  genet_grid_count(
    select_genet(
      genet_all = gene,
      Marker = m,
      level = tolower(l),
      Name = n,
      filters = filters
    )
  )
}

genet_grid_read <- function(config, filters){
  # Create grid object
  grid <- list()
  
  # Attacht metadata
  grid$config <- config
  grid$filters <- filters
  
  # Build it dataframe from metadata
  # From csv condition
  if(config$from_csv){
    grid$df <- read_genet_grid_csv(
      config$path, config$Marker, config$Level, config$Name
    )
  } else {
    dna <- load_edna_data(config$rarefied)
    
    grid$df <- build_genet_grid(
      config$Marker, config$Level, config$Name, filters
    )
  }
  
  # Total included
  if(!config$total){
    grid$df <- grid$df %>%
      filter(sample_loc != "TOTAL")
  }
  
  return(grid)
}




# ==== GRID READ ====

grid_read <- function(config, filters){
  
  # Initialize Grid
  grid <- init_grid(config, filters)
  
  # Fill Grid
  if(config$Type == "bota"){
    grid <- bota_grid_read(grid)
  }
  else if (config$Type == "genet"){
    grid <- genet_grid_read(grid)
  }
  
  return(grid)
}


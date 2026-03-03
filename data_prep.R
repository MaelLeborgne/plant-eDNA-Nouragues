
# ===== data_prep.R =====
# Regroup data preparation functions by Subject


# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)   # Basic
library(tibble)  # Tibble class
library(stringr) # Work on string
library(purrr)   # pmap_int()
library(vegan)   # OTU diversity

# ===== GENERAL =====

# Initialize empty grid
ref_grid <- function(){
  # X and Y
  grid <- expand.grid(
    x = seq(10, 290, 10), # X position of samples
    y = seq(10, 390, 10)  # Y position of samples
  )
  # Add indices i and j
  grid$i <- grid$x / 10
  grid$j <- grid$y / 10
  
  return(grid)
}

# Filtering the Data with the Name for a given Level
select_tax <- function(Data, Level, Name){
  
  # Select the focus in Data
  focus <- Data[Data[Level] == Name,]
  
  # Return
  return(focus)
}

# ===== BOTA =====

# Normalisation of column count of a grid
normalization <- function(count) {
  # Normalisation
  density <-  count / sum(count)
  # Return
  return(density)
}

bota_grid_engine <- function(bota_focus, 
                             R = 20,
                             fun_pixel,
                             var_name = "value",
                             normalize = TRUE){
  # Ref grid
  grid <- ref_grid()
  # Initailisation
  grid[[var_name]] <- NA
  
  # Iterate over pixels
  for (pixel in 1:nrow(grid)) {
    
    px <- grid$x[pixel]
    py <- grid$y[pixel]
    
    dij <- sqrt((px - bota_focus$Xfield)^2 + 
                  (py - bota_focus$Yfield)^2)
    
    sel <- which(dij < R & dij > 0)
    
    if (length(sel) == 0) {
      grid[[var_name]][pixel] <- 0
    } else {
      grid[[var_name]][pixel] <- fun_pixel(bota_focus, sel, dij)
    }
  }
  
  if (normalize) {
    grid[[paste0(var_name, "_density")]] <- 
      normalization(grid[[var_name]])
  }
  
  return(grid)
}

bota_grid_count <- function(bota_focus, R = 10){
  
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Count",
    fun_pixel = function(df, sel, dij){
      length(sel)
    }
  )
}

bota_grid_NCI <- function(bota_focus, R = 20){
  
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "NCI",
    fun_pixel = function(df, sel, dij){
      sum((df$CircCorr[sel]^2) * 
            exp(-0.25 * dij[sel]),
          na.rm = TRUE)
    }
  )
}

bota_grid_Gauss <- function(bota_focus, BETA, GAMMA, R = 20){
  
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Gaussian",
    fun_pixel = function(df, sel, dij){
      sum(
        (pi*(df$CircCorr[sel]/2)^2) *
          exp(-BETA * dij[sel]^GAMMA),
        na.rm = TRUE
      )
    }
  )
}

### Fonction to minimize (gauss)
rss_function <- function(par, ADNe_obs) {
  # Prameters
  beta <- par[1]
  gamma <- par[2]

  # Prediction
  pred_grid <- bota_grid_Gauss(beta, gamma)

  # Residuals
  resid <- ADNe_obs - pred_grid$density

  # RSS
  print(sum(resid^2))
  print(beta)
  sum(resid^2)
  }

### Launch Optimization (gauss)
optim_gauss <- function(genet_grid) {
  res <- optim(
    par = c(0.2,1),
    fn = rss_function,
    ADNe_obs = genet_grid$density,
    method = "L-BFGS-B",
    lower = c(0,0),  # logical constraint (biological)
  )
  return(res)
}
# OTU
bota_grid_OTU <- function(bota_focus, R = 20){
  
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "OTU",
    normalize = FALSE,   # pas de densité ici
    fun_pixel = function(df, sel, dij){
      length(unique(df$Species[sel]))
    }
  )
}
# Div_alpha
bota_grid_div_a <- function(bota_focus, R = 20){
  
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Div_alpha",
    normalize = FALSE,
    fun_pixel = function(df, sel, dij){
      sp_counts <- table(df$Species[sel])
      vegan::diversity(sp_counts, index = "simpson")
    }
  )
}


# Apply conditions on bota grid computation
bota_grid_compute <- function(bota, level, name, var_bota){
  
  # Sélection du taxon et gestion du all
  df <- if(name == "All") bota else select_tax(bota, level, name)
  
  # Retourne une grille vide si aucun individu
  if(nrow(df) == 0){
    warning("No individuals for this selection.")
    return(ref_grid())
  }
  
  # Choix de la métrique botanique via switch
  df <- switch(var_bota,
               "Count"     = bota_grid_count(df, R = 10),
               "NCI"       = bota_grid_NCI(df),
               "Gaussian"  = {
                 if(is.null(df)) stop("bota required for Gaussian")
                 result <- optim_gauss(bota)
                 bota_grid_Gauss(df, BETA = result$par[1], GAMMA = result$par[2])
               },
               "OTU"       = bota_grid_OTU(df, R = 20),
               "Div_alpha" = bota_grid_div_a(df, R = 20),
               stop("Unknown var_bota: ", var_bota)  # gestion d'erreur
  )
  
  return(df)
}

# ===== GENET =====

# Get Familly, Genus, species from faxid reference
get_FGs <-function(data){
  data <- left_join(data, tax_table, by="taxid")
  return(data)
}

# Select genetical reads with Marker and Name at a give Level
select_genet <- function(gene, Marker, level, Name){
  
  genet_focus <- gene[[Marker]] %>%          # Select marker
    get_FGs() %>%                            # Get Familly, Gender, species cols
    select_tax(level, Name) %>% # Select taxa
    filter(!is.na(.data[[level]]))           # Remove NA
  
  # Return
  return(genet_focus)
}

# Create genet grid for quantitative try
genet_grid_count <- function(genet_focus){
  
  # Create Grid
  grid <- ref_grid()
  
  # Get count of sequences by samples
  sample_df <- sample_seq_count(genet_focus)
  
  # Mean sample count by localisation
  genet_grid <- loc_mean(sample_df)
  
  # Normalisation
  genet_grid$Weight_density <- genet_grid$Weight %>% normalization()
  
  # Return
  return(genet_grid)
}

genet_grid_write <- function(PROJECT_PATH, focus, fileName){
  # Make grid
  grid <- genet_grid_count(focus)
  
  # Write
  write.csv(grid, paste0(PROJECT_PATH,"genet_grid/", fileName))
}

give_level_names <- function(Marker, Level){
  level <- tolower(Level)
  # Focus marker
  mark_focus <- gene[[Marker]] %>% # Select the marker
    get_FGs() %>% # Get Familly, Gender, species cols
    filter(!is.na(.data[[level]])) # remove NA
  # names in the level
  return(unique(mark_focus[[level]]))
}

genet_grid_write_all <- function(PROJECT_PATH, gene, markers, m_levels){
  # For each markers
  for (m in markers) {
    # Focus marker
    mark_focus <- gene[[m]] %>% # Select the marker
      get_FGs() # Get Familly, Gender, species cols
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
        genet_grid_write(PROJECT_PATH, focus, fileName)
      }
    }
  }
}

# Read grid from csv in ./genet_grid/
genet_grid_read <- function(PROJECT_PATH, m, l, n){
  fileName <- paste0(m,"_",l,"_",n,".csv")
  genet_grid <- read.csv(paste0(PROJECT_PATH,"genet_grid/",fileName))
  # Return
  return(genet_grid)
}

# ===== COMPARE =====

# Join genet and bota grid in the same dataframe
df_compare_join <- function(genet_grid, bota_grid, var_bota, var_genet){
  # Joining
  df_compare <- genet_grid %>%
    rename(value_genet = !!var_genet) %>%   # rename chosen genetic variable
    left_join(bota_grid, by = c("x", "y")) %>%
    rename(value_bota = !!var_bota)         # rename chosen botanical variable
  # Add error column
  df_compare <- df_compare %>%
    mutate(err_abs = value_genet - value_bota)
  # Return
  return(df_compare)
}

# Calcul de R^2 manuellement
R_squared <- function(df_compare, var_obs, var_pred){
  
  y <- df_compare[[var_obs]]                       # Donnée de référence
  y_pred <- df_compare[[var_pred]]                 # Prédictions du modèle
  
  y_mean <- mean(y, na.rm = TRUE)                  # Moyenne des données observées
  
  SS_residual <- sum((y - y_pred)^2, na.rm = TRUE) # Somme des carrés des résidus
  SS_total <- sum((y - y_mean)^2, na.rm = TRUE)    # Somme des carrés totaux
  
  R2 <- 1 - (SS_residual / SS_total)
  
  return(R2)
}

# Return Tax dataframe with all R2
R2_taxon <- function(Marker, Level, var_bota, var_genet){
  # All Names
  all_names <- give_level_names(Marker,Level)
  # Give all R2
  all_r2 <- lapply(all_names, function(n){
    # Bota grid
    bota_df <- bota_grid_compute(
      bota = bota,
      level = Level,
      name = n,
      var_bota = var_bota,
      var = var_genet
    )
    # Genet grid
    genet_df <- genet_grid_read(
      PROJECT_PATH,
      m = Marker,
      l = Level,
      n = n,
      var = var_genet
    )
    
    # Compute R2
    r2 <- df_compare_join(genet_df, bota_df, var_bota, var_genet) %>%
      R_squared(
        var_obs = var_genet, 
        var_pred = var_bota
      )
    
    # Put in a dataframe
    mono_df <- data.frame(
      Name = n, # name
      R2 = round(r2,digits = 2)
    )
    return(mono_df)
  })
  
  # Fusion en un seul dataframe
  final_df <- do.call(rbind, all_r2) %>% 
    filter(!is.na(.data[["R2"]])) %>%
    arrange(desc(R2))
  
  return(final_df)
}

# ==== WRITE ALL GENET_GRID ====
# # Test loop
# markers <- c("GH","18S","ITS")
# m_levels <- c("Family","Genus","Species")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)

# # Test write (Count)
# markers <- c("GH")
# m_levels <- c("Family")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels, "Count")

# Test write (Div_alpha)
# markers <- c("GH")
# m_levels <- c("Family")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)



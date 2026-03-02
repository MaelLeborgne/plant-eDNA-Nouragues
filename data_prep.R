
# ===== data_prep.R =====
# Regroup data preparation functions by Subject


# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # Basic
library(tibble) # Tibble class
library(stringr) # Work on string
library(purrr) # pmap_int()

# ===== GENERAL =====

# Initialize empty grid
ref_grid <- function(){
  # X and Y
  grid <- expand.grid(
    x = seq(10, 300, 10), # X position of samples
    y = seq(10, 400, 10)  # Y position of samples
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

# Make bota grid (count and density)
bota_grid_count <- function(bota_focus, R){
  
  # Create Grid
  grid <- ref_grid()
  
  # Initialize Count to NA
  grid$count <- NA
  
  # Count the number of trees around each pixels in a circle of radius R
  grid$count <- pmap_int( # Apply a function to each coordinates (x,y)
    list(grid$x, grid$y), # Element of iterations
    \(x, y) sum((x - bota_focus$Xfield)^2 + # Count TRUE ...
                  (y - bota_focus$Yfield)^2 <= R^2, # ... on Pytagore condition (dx²-dy² <= R²)
                na.rm = TRUE)
  )
  
  # Convert count to density
  grid$density <- grid$count %>% normalization()
  
  # Return
  return(grid)
}

bota_grid_NCI <- function(bota_focus){
  
  # Create Grid
  grid <- ref_grid()
  
  # Initialize NCI to NA
  grid$count <- NA
  
  # Loop over rows (pixels)
  for (pixel in 1:nrow(grid)) {
    
    # Coordinates of the current pixel
    px <- grid$x[pixel]
    py <- grid$y[pixel]
    
    # Euclidean distance between this pixel and all trees
    dij <- sqrt((px - bota_focus$Xfield)^2 + 
                  (py - bota_focus$Yfield)^2)
    
    # Keep only trees within 20 m (excluding distance 0)
    sel <- which(dij < 20 & dij > 0)
    
    if (length(sel) == 0) {
      grid$count[pixel] <- 0
    } else { # Compute the NCI
      grid$count[pixel] <- sum( 
        (bota_focus$CircCorr[sel]^2) * 
          exp(-0.25 * dij[sel]),
        na.rm = TRUE
      )
    }
  }
  
  # Convert NCI to density
  grid$density <- grid$count %>% normalization()
  
  # Return
  return(grid)
}
  
# Gaussian model for botanical grid
bota_grid_Gauss <- function(bota_focus, BETA, GAMMA){
  
  # Create Grid
  grid <- ref_grid()
  # Initialize NCI to NA
  grid$count <- NA
  
  # Loop over rows (pixels)
  for (pixel in 1:nrow(grid)) {
    
    # Coordinates of the current pixel
    px <- grid$x[pixel]
    py <- grid$y[pixel]
    
    # Euclidean distance between this pixel and all trees
    dij <- sqrt((px - bota_focus$Xfield)^2 + 
                  (py - bota_focus$Yfield)^2)
    
    # Keep only trees within 20 m (excluding distance 0)
    sel <- which(dij < 20 & dij > 0)
    
    if (length(sel) == 0) {
      grid$count[pixel] <- 0
    } else { # Compute the gaussian index
      grid$count[pixel] <- sum( 
        (pi*(bota_focus$CircCorr[sel]/2)^2) * # Basal Area
          exp(-BETA * dij[sel]^GAMMA),              # Gaussian Kernel
        na.rm = TRUE
      )
    }
  }
  
  # Convert NCI to density
  grid$density <- grid$count %>% normalization()
  
  # Return
  return(grid)
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

### Launch Optimization
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

# Create genet grid
genet_grid_count <- function(genet_focus){
  
  # Create Grid
  grid <- ref_grid()
  
  # Get count of sequences by samples
  sample_df <- sample_seq_count(genet_focus)
  
  # Mean sample count by localisation
  genet_grid <- loc_seq_count(sample_df)
  
  # Normalisation
  genet_grid$density <- genet_grid$count %>% normalization()
  
  # Return
  return(genet_grid)
}

genet_grid_divaplha <- function(genet_focus){
  # Create Grid
  grid <- ref_grid()
  
  # Get diversity of sequences by localisation
  genet_grid <- ???
  
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
df_compare_join <- function(genet_grid, bota_grid){
  # Joining
  df_compare <- genet_grid %>%
    rename(dens_genet = density) %>% # rename
    rename(count_genet = count) %>% # rename
    left_join(bota_grid, by = c("x", "y")) %>% # left join
    rename(dens_bota = density) %>% # rename density
    rename(count_bota = count) # rename
  # Add error column
  df_compare <- df_compare %>%
    mutate(err_abs = dens_genet - dens_bota)
  # Return
  return(df_compare)
}

R_squared <- function(df_compare){
  # Calcul de R^2 manuellement
  y <- df_compare$dens_genet             # Donnée de référence
  y_pred <- df_compare$dens_bota         # Prédictions du modèle
  y_mean <- mean(df_compare$dens_genet)  # Moyenne des données observées
  SS_residual <- sum((y - y_pred)^2)  # Somme des carrés des résidus
  SS_total <- sum((y - y_mean)^2)    # Somme des carrés totaux
  
  # Return
  R2 <- 1 - (SS_residual / SS_total)
  return(R2)
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



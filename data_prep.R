
# ===== data_prep.R =====
# Regroup data preparation functions by Subject


# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # Basic
library(tibble) # Tibble class
library(stringr) # Work on string

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
genet_grid <- function(genet_focus){
  
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

genet_grid_write <- function(PROJECT_PATH, focus, fileName){
  # Make grid
  grid <- genet_grid(focus)
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
      l_names <- unique(mark_focus[[l]])
      # For each names
      for (n in l_names) {
        N <- str_to_sentence(n) # upper first letter
        # File name
        fileName <- paste0(m,"_",L,"_",N,".csv")
        print(fileName)
        # Focus taxa
        focus <- select_tax(mark_focus, l, n) %>%
          filter(!is.na(.data[[l]])) # col family, genus or specie
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

# ==== WRITE ALL GENET_GRID ====
# # Test loop
# markers <- c("GH","18S","ITS")
# m_levels <- c("Family","Genus","Species")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)

# Test write
# markers <- c("GH")
# m_levels <- c("Family")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)



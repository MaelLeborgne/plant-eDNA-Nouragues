
# ===== data_prep.R =====
# Regroup data preparation functions by Subject

# Path of the project
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"

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

# Normalisation of column count of a grid
normalization <- function(count) {
  # Normalisation
  density <-  count / sum(count)
  # Return
  return(density)
}

# Generic function to compute bota grid
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
    grid[[var_name]] <- normalization(grid[[var_name]])
  }
  
  return(grid)
}

# Count
bota_grid_count <- function(bota_focus = bota, 
                            R = 10, 
                            user_norm = FALSE
                            ){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Count",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      length(sel)
    }
  )
}

# NCI
bota_grid_NCI <- function(bota_focus, R = 20, user_norm){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "NCI",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      sum((df$CircCorr[sel]^2) * 
            exp(-0.25 * dij[sel]),
          na.rm = TRUE)
    }
  )
}

# Gaussian
bota_grid_Gauss <- function(bota_focus, BETA, GAMMA, R = 20, user_norm){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Gaussian",
    normalize = user_norm,
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
bota_grid_OTU <- function(bota_focus, R = 20, user_norm){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "OTU",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      length(unique(df$Species[sel]))
    }
  )
}

# Div_alpha
bota_grid_div_a <- function(bota_focus, R = 20, user_norm){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Div_alpha",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      sp_counts <- table(df$Species[sel])
      vegan::diversity(sp_counts, index = "simpson")
    }
  )
}

# Apply conditions for bota grid calcul
bota_grid_compute <- function(bota_df   = bota, 
                              level     = "Family", 
                              name      = "Fabaceae", 
                              var_bota  = "Count", 
                              user_norm = TRUE,
                              radius = 15
                              ){
  
  # Sélection du taxon et gestion du all
  df <- if(name == "All") bota_df else select_tax(bota_df, level, name)
  
  # # Retourne une grille vide si aucun individu
  # if(nrow(df) == 0){
  #   warning("No individuals for this selection.")
  # }
  
  # Choix de la métrique botanique via switch
  df <- switch(var_bota,
               "Count"     = bota_grid_count(df, R = radius, user_norm),
               "NCI"       = bota_grid_NCI(df, user_norm),
               "Gaussian"  = {
                 if(is.null(df)) stop("bota required for Gaussian")
                 result <- optim_gauss(bota)
                 bota_grid_Gauss(df, BETA = result$par[1], GAMMA = result$par[2], user_norm)
               },
               "OTU"       = bota_grid_OTU(df, R = radius, user_norm),
               "Div_alpha" = bota_grid_div_a(df, R = radius, user_norm),
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
select_genet <- function(genet_all = gene, 
                         Marker = "GH", 
                         level = "family", 
                         Name = "All"){
  
  genet_focus <- gene[[Marker]] %>%          # Select marker
    get_FGs() %>%                            # Get Familly, Gender, species cols
    select_tax(level, Name) %>% # Select taxa
    filter(!is.na(.data[[level]]))           # Remove NA
  
  # Return
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

genet_grid_write_all <- function(path = PROJECT_PATH, 
                                 genet_all = gene, 
                                 markers  ="GH", 
                                 m_levels = "Family"){
  # For each markers
  for (m in markers) {
    # Focus marker
    mark_focus <- genet_all[[m]] %>% # Select the marker
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
        genet_grid_write(path, focus, fileName)
      }
    }
  }
}

# Read grid from csv in ./genet_grid/
genet_grid_read <- function(path = PROJECT_PATH,
                            m = "GH", 
                            l = "Family", 
                            n = "Fabaceae"
                            ){
  fileName <- paste0(m,"_",l,"_",n,".csv")
  genet_grid <- read.csv(paste0(path,"genet_grid/",fileName))
  # Return
  return(genet_grid)
}

# ===== COMPARE =====

# Join genet and bota grid in the same dataframe
df_compare_join <- function(grid_pred  = genet_grid_read(), 
                            grid_obs   = bota_grid_compute(), 
                            var_pred   = "Weight_density", 
                            var_obs    = "Count"
                            ){
  # Joining
  df_compare <- grid_pred %>%
    rename(value_pred = !!var_pred) %>%   # rename chosen genetic variable
    left_join(grid_obs, by = c("x", "y")) %>%
    rename(value_obs = !!var_obs)         # rename chosen botanical variable
  # Add error column
  df_compare <- df_compare %>%
    mutate(err_abs = value_pred - value_obs)
  # Return
  return(df_compare)
}

# Calcul de R^2 manuellement
R_squared <- function(df_compare){
  
  y <- df_compare[["value_obs"]]                       # Donnée de référence
  y_pred <- df_compare[["value_pred"]]                 # Prédictions du modèle
  
  y_mean <- mean(y, na.rm = TRUE)                  # Moyenne des données observées
  
  SS_residual <- sum((y - y_pred)^2, na.rm = TRUE) # Somme des carrés des résidus
  SS_total <- sum((y - y_mean)^2, na.rm = TRUE)    # Somme des carrés totaux
  
  R2 <- 1 - (SS_residual / SS_total)
  
  return(R2)
}

# Return Tax dataframe with all R2
R2_taxon <- function(Marker    = "GH", 
                     Level     = "Family", 
                     var_bota  = "Count", 
                     var_genet = "Weight_density", 
                     user_norm = TRUE
                     ){
  # All Names
  all_names <- give_level_names(Marker,Level)
  # Give all R2
  all_r2 <- lapply(all_names, function(n){
    # Bota grid
    bota_df <- bota_grid_compute(
      bota_df      = bota,
      level        = Level,
      name         = n,
      var_bota     = var_bota,
      user_norm    = user_norm
    )
    # Genet grid
    genet_df <- genet_grid_read(
      PROJECT_PATH,
      m = Marker,
      l = Level,
      n = n
    )
    
    # Compute R2
    r2 <- df_compare_join(
      grid_pred  = genet_df, 
      grid_obs   = bota_df, 
      var_pred   = var_genet, 
      var_obs    = var_bota
    ) %>% R_squared()
    
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

# # Test write
# markers <- c("GH")
# m_levels <- c("Family")
# genet_grid_write_all(PROJECT_PATH, gene, markers, m_levels)



# ==== RICHNESS TABLE ====

summary_per_clade <- function(Marker = "GH",
                              Level  = "Family",
                              fun_samples = function(df) df %>% filter(sample_loc != "TOTAL") %>% pull(OTU),
                              fun_total   = function(df) df %>% filter(sample_loc == "TOTAL") %>% pull(OTU),
                              digits = 0,
                              per_sample = TRUE) {
  
  # Récupère tous les noms du clade
  all_names <- give_level_names(Marker, Level)
  # Grid with all names
  genet_df_all <- genet_grid_read(n="All")
  
  results <- lapply(all_names, function(n) {
    
    # Lit les données
    genet_df <- genet_grid_read(PROJECT_PATH,
                                m = Marker,
                                l = Level,
                                n = n)
    
    # Extraie le vecteur selon la fonction fournie
    vec <- fun_samples(genet_df)
    
    # Total, moyenne, sd
    total <- fun_total(genet_df, genet_df_all)
    mean_val <- mean(vec, na.rm = TRUE)
    sd_val   <- sd(vec, na.rm = TRUE)
    
    # Formatage publication
    total_fmt <- format(round(total, digits),
                        big.mark = ",",
                        scientific = FALSE)
    mean_fmt <- round(mean_val, digits)
    sd_fmt   <- round(sd_val, digits)
    
    if(per_sample){
      data.frame(
        Clade = n,
        Metric = paste0(total_fmt, " (", mean_fmt, " ± ", sd_fmt, ")")
      )
    }
    else{
      data.frame(
        Clade = n,
        Metric = paste0(total_fmt)
      )
    }
    
  })
  
  final_df <- do.call(rbind, results)
  return(final_df)
}

# Join metrics
tab_OTU_stat <- function(Marker   = "GH",
                         Level    = "Family") {
  
  # Richness
  df_rich <- summary_per_clade(
    Marker = Marker,
    Level  = Level,
    fun_samples = function(df) df %>% filter(sample_loc != "TOTAL") %>% pull(OTU),
    fun_total   = function(df,df_t) df %>% filter(sample_loc == "TOTAL") %>% pull(OTU),
    digits = 0,
    per_sample = TRUE
  ) %>% rename(Richness = Metric)
  
  # Diversity
  df_div <- summary_per_clade(
    Marker = Marker,
    Level  = Level,
    fun_samples = function(df) df %>% filter(sample_loc != "TOTAL") %>% pull(Div_alpha),
    fun_total   = function(df, df_t) df %>% filter(sample_loc == "TOTAL") %>% pull(Div_alpha),
    digits = 1,
    per_sample = TRUE
  ) %>% rename(Diversity = Metric)
  
  # Relative abundance
  df_rel <- summary_per_clade(
    Marker = Marker,
    Level  = Level,
    fun_samples = function(df) NULL,
    fun_total = function(df, df_t) {
      W_clade <- df %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
      W_all <- df_t %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
      return(W_clade * 100 / W_all)
      },
    digits = 2,
    per_sample = FALSE
  ) %>% rename(RelAb = Metric)
  
  # Fusion
  df_join <- df_rich %>%
    left_join(df_div, by = "Clade") %>%
    left_join(df_rel, by = "Clade") %>%
    arrange(Clade)
  
  return(df_join)
}

# ==== OUT OF APP ====

plot_abrel_comp <- function(){
  library(ggplot2)
  # Genet
  gene_stat <- tab_OTU_stat()
  gene_stat$RelAb <- gene_stat$RelAb %>% as.numeric()
  # Bota
  tab <- table(bota$Family)
  df <- as.data.frame(tab)
  colnames(df) <- c("Clade", "Richesse")
  df$Abond_rel <- df$Richesse*100 / sum(df$Richesse)
  # Join
  df_rich_comp_f <- left_join(gene_stat, df, by="Clade") %>% 
    filter(RelAb >1) %>% select(Clade, RelAb, Abond_rel)
  # Plot
  ggplot(df_rich_comp_f,aes(x=RelAb,y=Abond_rel)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey70")
}

plot_abrel_BA <- function(){
  library(ggplot2)
  # Genet
  gene_stat <- tab_OTU_stat()
  gene_stat$RelAb <- gene_stat$RelAb %>% as.numeric()
  # Bota
  sum_BA <- bota %>% 
    mutate(BA_indiv = (bota$CircCorr^2) / (4 * pi * 10000)) %>% 
    group_by(Family) %>% 
    summarise(BA=sum(BA_indiv)) %>%
    rename(Clade = Family)
  # Join
  df_rich_comp_f <- left_join(gene_stat, sum_BA, by="Clade") %>% 
    filter(RelAb >1) %>% select(Clade, RelAb, BA)
  # Plot
  ggplot(df_rich_comp_f,aes(x=RelAb,y=BA)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey70")
}

plot_abrel_BA_genus <- function(){
  library(ggplot2)
  # Genet
  gene_stat <- tab_OTU_stat(Level = "Genus")
  gene_stat$RelAb <- gene_stat$RelAb %>% as.numeric()
  # Bota
  sum_BA <- bota %>% 
    mutate(BA_indiv = (bota$CircCorr^2) / (4 * pi * 10000)) %>% 
    group_by(Genus) %>% 
    summarise(BA=sum(BA_indiv)) %>%
    rename(Clade = Genus)
  # Join
  df_rich_comp_f <- left_join(gene_stat, sum_BA, by="Clade") %>% 
    filter(RelAb >0.01) %>% select(Clade, RelAb, BA)
  # Plot
  ggplot(df_rich_comp_f,aes(x=RelAb,y=BA)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey70")
}

plot_borras_2A <- function(){
  # Genet
  gene_df <- genet_grid_read(n = "All") %>% 
    filter(sample_loc != "TOTAL")
  # Bota
  bota_rad <- lapply(seq(5,100,5), function(R){
    bota_df <- bota_grid_compute(
      bota_df    = bota,
      name       = "All",
      var_bota   = "OTU",
      user_norm  = FALSE,
      radius = R
    ) %>% select()
    return(bota_df)
  })
  

  # Fusion en un seul dataframe
  final_df <- do.call(rbind, bota_rad) %>% 
    filter(!is.na(.data[["R2"]])) %>%
    arrange(desc(R2))
  
  return(final_df)
  
  
  
}
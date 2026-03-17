
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
library(ggplot2) # Plot


# ===== COMPARE =====

# Join genet and bota grid in the same dataframe
df_compare_join <- function(grid_pred  = genet_grid_read(), 
                            grid_obs   = bota_grid_read(), 
                            var_pred   = "Weight_density", 
                            var_obs    = "Count",
                            total = FALSE
                            ){
  # Joining
  df_compare <- grid_pred %>%
    rename(value_pred = !!var_pred) %>%   # rename chosen genetic variable
    left_join(grid_obs, by = c("x", "y")) %>%
    rename(value_obs = !!var_obs)         # rename chosen botanical variable
  # Add error column
  df_compare <- df_compare %>%
    mutate(
      err_abs = value_pred - value_obs
      )
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

corr_b_g <- function(df_compare){
  cor.test(df_compare[["value_obs"]], df_compare[["value_pred"]], method = "pearson")$estimate
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
    bota_df <- bota_grid_read(
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
      n = n,
      total = FALSE
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


# ==== RICHNESS TABLE ====

summary_per_clade <- function(Marker = "GH",
                              Level  = "Family",
                              fun_samples = function(df) df %>% pull(OTU),
                              fun_total   = function(df, df_t) {
                                W_clade <- df %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
                                W_all <- df_t %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
                                return(as.numeric((W_clade * 100) / W_all))
                              },
                              digits = 0,
                              per_sample = FALSE) {
  
  # Récupère tous les noms du clade
  all_names <- give_level_names(Marker, Level, filters = filter_list)
  # Grid with all names
  genet_df_all <- genet_grid_read(m = Marker, 
                                  n = "All",
                                  from_csv = TRUE,
                                  total = TRUE)
  
  results <- lapply(all_names, function(n) {
    
    # Lit les données
    genet_df <- genet_grid_read(PROJECT_PATH,
                                m = Marker,
                                l = Level,
                                n = n,
                                from_csv = TRUE,
                                total = TRUE)
    
    # Remove TOTAL
    genet_df_clean <- genet_df %>% filter(sample_loc != "TOTAL")
    
    # Extraie le vecteur selon la fonction fournie
    vec <- fun_samples(genet_df_clean)
    
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
        Metric = total
      )
    }
    
  })
  
  final_df <- do.call(rbind, results)
  return(final_df)
}

tab_rel_ab <- function(Marker   = "GH",
                       Level    = "Family"){
  
  # Compute table
  df_rel <- summary_per_clade(
    Marker = Marker,
    Level  = Level,
    fun_samples = function(df) NULL,
    fun_total = function(df, df_t) {
      W_clade <- df %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
      W_all <- df_t %>% filter(sample_loc == "TOTAL") %>% pull(Weight)
      return(as.numeric((W_clade * 100) / W_all))
    },
    digits = 2,
    per_sample = FALSE
  ) %>% rename(Genet_ab_rel = Metric)
  
  return(df_rel)
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
    fun_samples = function(df) df %>% filter(sample_loc != "TOTAL") %>% pull(Shannon),
    fun_total   = function(df, df_t) df %>% filter(sample_loc == "TOTAL") %>% pull(Shannon),
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
      return(as.numeric((W_clade * 100) / W_all))
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

# ==== TRY : GROWTH_FORM ====
# Commende linux : awk -F'\t' '$10=="42" {print $5 "\t" $15}' 47883.txt | uniq > growth_forms.txt

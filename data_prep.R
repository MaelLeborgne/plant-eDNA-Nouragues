
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
  
  grid <- left_join(grid, samples_ref, by = c("x","y"))
  
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
bota_grid_engine <- function(bota_focus = bota, 
                             R = 20,
                             fun_pixel = function(df, sel, dij){length(sel)},
                             var_name = "Count",
                             normalize = TRUE,
                             return_comm = FALSE){
  
  grid <- ref_grid()
  
  if(return_comm){
    grid[[var_name]] <- vector("list", nrow(grid))
  } else {
    grid[[var_name]] <- NA_real_
  }
  
  for (pixel in 1:nrow(grid)) {
    # adresse du pixel
    px <- grid$x[pixel]
    py <- grid$y[pixel]
    # distance de tout les arbre avec le pixel
    dij <- sqrt((px - bota_focus$Xfield)^2 + 
                  (py - bota_focus$Yfield)^2)
    # Sélectionner les arbres dans un intervale de distance
    sel <- which(dij < R & dij > 0)
    
    if (length(sel) == 0) {
      
      if(return_comm){
        grid[[var_name]][[pixel]] <- table()
      } else {
        grid[[var_name]][pixel] <- 0
      }
      
    } else {
      
      if(return_comm){
        grid[[var_name]][[pixel]] <- fun_pixel(bota_focus, sel, dij)
      } else {
        grid[[var_name]][pixel] <- fun_pixel(bota_focus, sel, dij)
      }
      
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

grid_basal_area <- function(bota_focus = bota, 
                            R = 10, 
                            user_norm = FALSE
){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Basal_area_m2",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      sum((bota[sel, ]$CircCorr^2) / (4 * pi * 10000))
    }
  )
}

# NCI
bota_grid_NCI <- function(bota_focus, R = 20, user_norm =FALSE){
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
bota_grid_Gauss <- function(bota_focus, BETA, GAMMA, R = 20, user_norm=FALSE){
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
bota_grid_OTU <- function(bota_focus, R = 15, user_norm){
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

bota_grid_OTU_matrice <- function(bota_focus, R = 15){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "OTU",
    normalize = FALSE,
    return_comm = TRUE, # Matrice
    fun_pixel = function(df, sel, dij){
      table(df$Species[sel])
    }
  )
}

# Shannon
bota_grid_div_a <- function(bota_focus, R = 20, user_norm){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "Shannon",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      sp_counts <- table(df$Species[sel])
      vegan::diversity(sp_counts, index = "simpson")
    }
  )
}

# Shannon
bota_grid_H1 <- function(bota_focus, R = 20, user_norm=FALSE){
  # Give calcul to generic function
  bota_grid_engine(
    bota_focus,
    R = R,
    var_name = "H1",
    normalize = user_norm,
    fun_pixel = function(df, sel, dij){
      sp_counts <- table(df$Species[sel])
      shannon <- vegan::diversity(sp_counts, index = "simpson")
      exp(shannon)
    }
  )
}

# Apply conditions for bota grid calcul
bota_grid_compute <- function(bota_df   = bota, 
                              level     = "Family", 
                              name      = "All", 
                              var_bota  = "Count", 
                              user_norm = TRUE,
                              return_comm = FALSE,
                              radius = 15
                              ){
  
  # Sélection du taxon et gestion du all
  df <- if(name == "All") bota_df else select_tax(bota_df, level, name)
  
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
               "OTU_mat"   = bota_grid_OTU_matrice(df, R = radius),
               "Shannon" = bota_grid_div_a(df, R = radius, user_norm),
               stop("Unknown var_bota: ", var_bota)  # gestion d'erreur
  )
  
  return(df)
}



# ===== GENET =====

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

genet_grid_write <- function(PROJECT_PATH, focus, fileName){
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

genet_grid_write_all <- function(path = PROJECT_PATH, 
                                 genet_all = gene, 
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
        genet_grid_write(path, focus, fileName)
      }
    }
  }
}

# Read grid from csv in ./genet_grid/
genet_grid_read <- function(path = PROJECT_PATH,
                            m = "GH", 
                            l = "Family", 
                            n = "Fabaceae",
                            from_csv = TRUE,
                            rarefied = FALSE,
                            filters = list(),
                            total = FALSE
                            ){
  
  # Init
  genet_grid <- data.frame()
  
  if(from_csv){ # Read from csv
    fileName <- paste0(m,"_",l,"_",n,".csv")
    genet_grid <- read.csv(paste0(path,"genet_grid_5filters/",fileName))
  }
  else{ # Read from
    if(rarefied){
      # ---- eDNA ----
      dna_comm <- readRDS("/home/jean/Bureau/GH_clean/GH_mblist_annot_metabar_clean_clust97_RAR.rds")
      dna_comm <- dna_comm$reads
      dna_rich <- rowSums(dna_comm > 0)
    }
    genet_grid <- genet_grid_count(
      select_genet(
        genet_all = gene, 
        Marker = m, 
        level = tolower(l), 
        Name = n,
        filters = filters)
    )
  }
  
 
  if (!total){ # remove de total
    genet_grid <- genet_grid   %>% 
      filter(sample_loc != "TOTAL")
  }

  # Return
  return(genet_grid)
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

# ==== TRY : GROWTH_FORM ====
# Commende linux : awk -F'\t' '$10=="42" {print $5 "\t" $15}' 47883.txt | uniq > growth_forms.txt

# ===== COMPARE =====

# Join genet and bota grid in the same dataframe
df_compare_join <- function(grid_pred  = genet_grid_read(), 
                            grid_obs   = bota_grid_compute(), 
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
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  
  # ---- Genet ----
  gene_df <- genet_grid_read(m = "GH", 
                             n = "All", 
                             total = FALSE) %>% 
    select(x, y, OTU) %>%
    rename(DNA = OTU)
  
  # ---- Bota ----
  radii <- seq(5, 100, 5)
  
  bota_rad <- map(radii, function(R){
    
    bota_grid_compute(
      bota_df   = bota,
      name      = "All",
      var_bota  = "OTU",
      user_norm = FALSE,
      radius    = R
    ) %>% 
      select(x, y, OTU) %>%
      rename(value = OTU) %>%
      mutate(method = as.character(R))
  })
  
  bota_rad <- bind_rows(bota_rad)
  
  # ---- Join ----
  # Ajouter DNA
  gene_long <- gene_df %>%
    rename(value = DNA) %>%
    mutate(method = "DNA")
  
  # Fusion long
  final_long <- bind_rows(gene_long, bota_rad)
  
  # Garcer l'ordre
  final_long$method <- factor(
    final_long$method,
    levels = c("DNA", as.character(radii))
  )
  
  # ---- Boxplot ----
  ggplot(final_long, aes(x = method, y = value)) +
    geom_boxplot() +
    labs(x = "Radius / Method", y = "OTU richness") +
    theme_bw()
  
}

plot_borras_2E <- function(){ # TODO : change source of genet data

  library(dplyr)
  library(purrr)
  library(vegan)
  library(ggplot2)
  
  radii <- seq(10,100,5)
  
  cor_results <- map(radii, function(rad){
    
    # ---- Stem census (bota) ----
    bota_df <- bota_grid_engine(
      R = rad,
      var_name = "OTU",
      normalize = FALSE,
      return_comm = TRUE, # Matrice
      fun_pixel = function(df, sel, dij){
        if(length(sel) == 0) return(integer(0))
        table(df$Species[sel])
      }
    )
    
    stem_comm <- bind_rows(
      lapply(seq_len(nrow(bota_df)), function(i){
        
        tab <- bota_df$OTU[[i]]
        if(length(tab) == 0) return(NULL)
        
        df <- as.data.frame(t(as.numeric(tab)))
        names(df) <- names(tab)
        
        df$sample_id <- bota_df$sample_id[i]
        
        df
      })
    )
    
    stem_comm[is.na(stem_comm)] <- 0
    
    # mettre sample_id en rownames
    rownames(stem_comm) <- stem_comm$sample_id
    stem_comm$sample_id <- NULL
    
    # matrice
    stem_comm <- as.matrix(stem_comm)
    
    # richesse
    stem_rich <- rowSums(stem_comm > 0)
    
    min_n <- min(stem_rich[stem_rich > 0])
    keep <- stem_rich >= min_n
    stem_comm_filtered <- stem_comm[keep, , drop = FALSE]
    
    # rarefaction basée sur le nombre d'individus
    stem_rare <- vegan::rarefy(stem_comm_filtered, sample = min_n)
    
    # ---- eDNA ----
    dna_comm <- readRDS("/home/jean/Bureau/GH_clean/GH_mblist_annot_metabar_clean_clust97_RAR.rds")
    dna_comm <- dna_comm$reads
    dna_rich <- rowSums(dna_comm > 0)
    
    # ---- Intersect ----
    common <- intersect(names(stem_rich), names(dna_rich))
    
    # ---- Pearson correlation ----
    test <- cor.test(
      stem_rich[common],
      dna_rich[common],
      method = "pearson"
    )
    
    # return
    tibble(
      radius = rad,
      r      = test$estimate,
      CI_low = test$conf.int[1],
      CI_high= test$conf.int[2],
      p      = test$p.value
    )
    
  })
  
  # ---- Bind ----
  cor_results <- bind_rows(cor_results)
  
  # ---- Plot ----
  
  ggplot(cor_results, aes(x = radius, y = r)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(
      size = 4,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 1
    ) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      width = 1,
      colour = "black"
    ) +
    scale_x_continuous(
      breaks = seq(10, 100, 10)
    ) +
    theme_bw() +
    labs(
      x = "Radius (m)",
      y = "Pearson correlation",
      title = "Stem OTU vs DNA OTU richness across spatial scales"
    )
}

plot_bota_intro <- function(){
  library(ggplot2)
  
  # Supposons que ref_grid() renvoie un data.frame avec X et Y
  grid_df <- ref_grid()  # doit avoir colonnes X, Y
  
  ggplot(data = bota, aes(x = Xfield, y = Yfield)) +
    
    # points principaux
    geom_point(
      aes(color = .data[["Family"]], size = CircCorr),
      shape = 21,
      fill = NA,
      stroke = 0.5
    ) +
    
    # croix noires pour la grille
    geom_point(
      data = grid_df,
      aes(x = x, y = y),
      shape = 4,    # shape 4 = croix
      color = "black",
      size = 2
    ) +
    
    # cadre autour de la parcelle
    geom_rect(
      aes(xmin = 0, xmax = 300, ymin = 0, ymax = 400),
      fill = NA,
      color = "black",
      linewidth = 0.5
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(size = 16)
      ) +
    coord_fixed()
}

plot_Fabaceae_ex <- function(){
  plot_grid(genet_grid_read(m ="GH", l="Family", n="Fabaceae"), 
            var="Weight", 
            stems = TRUE, 
            bota_f = bota %>% filter(Family=="Fabaceae")) +
    theme(text = element_text(size = 16))
}

plot_tab_familly <- function(){
  
  # Genet
  tab_genet <- tab_OTU_stat()
  tab_genet$RelAb <- as.numeric(tab_genet$RelAb)
  # Other : calcul
  other_RelAb <- tab_genet %>% 
    select(RelAb) %>% 
    filter(RelAb < 1) %>% 
    sum()
  # Other : line
  genet_other <- data.frame(
    Clade = "Other",
    RelAb = other_RelAb
  )
  # Join >1% and other
  tab_genet <- tab_genet %>% 
    filter(RelAb >= 1) %>%
    bind_rows(genet_other) %>%
    select(Clade, RelAb) %>%
    rename(Genet_ab_rel = RelAb)
  
  # Bota
  tab_bota <- table(bota$Family)
  tab_bota <- as.data.frame(tab_bota)
  colnames(tab_bota) <- c("Clade", "Bota_ab_rel")
  tab_bota$Bota_ab_rel <- tab_bota$Bota_ab_rel*100 / sum(tab_bota$Bota_ab_rel)
  # Other : calcul
  other_RelAb <- tab_bota %>% 
    select(Bota_ab_rel) %>% 
    filter(Bota_ab_rel < 1) %>% 
    sum()
  # Other : line
  bota_other <- data.frame(
    Clade = "Other",
    Bota_ab_rel = other_RelAb
  )
  # Join >1% and other
  tab_bota <- tab_bota %>% 
    filter(Bota_ab_rel >= 1) %>%
    bind_rows(bota_other) %>%
    select(Clade, Bota_ab_rel)
  # Round
  tab_bota$Bota_ab_rel <- tab_bota$Bota_ab_rel
  
  # Join
  df_rel_ab <- inner_join(tab_genet, tab_bota, by="Clade") %>%
    arrange(ifelse(Clade == "Other", 1, 0), desc(Genet_ab_rel))
  
  # # Replace NA by 0
  # df_rel_ab <- df_rel_ab %>%
  #   mutate(Bota_ab_rel = replace_na(Bota_ab_rel, 0)) %>%
  #   mutate(Genet_ab_rel = replace_na(Genet_ab_rel, 0))
  
  # Mettre en format long
  df_long <- df_rel_ab %>%
    pivot_longer(
      cols = c(Genet_ab_rel, Bota_ab_rel),
      names_to = "Type",
      values_to = "Abondance"
    )
  
  # Spearman correlation with p-value
  spearman_test <- cor.test(df_rel_ab$Genet_ab_rel,
                            df_rel_ab$Bota_ab_rel,
                            method = "spearman")
  
  # Résultat
  rho <- spearman_test$estimate   # rho
  pval <- spearman_test$p.value    # p-value
  
  # Formater la p-value en texte lisible
  p_text <- if (pval < 0.001) {
    "< 0.001"
  } else if (pval < 0.01) {
    "< 0.01"
  } else if (pval < 0.05) {
    "< 0.05"
  } else if (pval < 0.1) {
    "< 0.1"
  } else {
    "n.s."
  }
  
  # Barplot horizontal avec texte de corrélation
  ggplot(df_long, aes(x = Abondance, y = reorder(Clade, Abondance), fill = Type)) +
    geom_col(position = "dodge") +
    labs(
      x = "Abondance relative (%)",
      y = "Clade",
      fill = "Type d'inventaire",
      title = "Richesse taxonomique relative par famille selon l'inventaire Botanique et Génétique (GH)"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Genet_ab_rel" = "#1f78b4", "Bota_ab_rel" = "#33a02c")) +
    annotate(
      "text",
      x = max(df_long$Abondance) * 0.7,  # position horizontale
      y = 5,                            # position verticale (au-dessus du bas du plot)
      label = paste0("Spearman 'rho = ", round(rho, 2), "\n p",p_text),
      hjust = 0.5,                          # aligné à gauche
      size = 3.5,
      color = "black"
    ) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.line = element_line(color = "gray"),  # Lignes des axes en noir
      axis.text = element_text(color = "black"),  # Texte des axes en noir
      axis.ticks = element_line(color = "black")  # Ticks des axes en noir
    )
  
}


# ggplot(df, aes(x = group, fill = category)) +
#   geom_bar(position = "fill")

plot_tab_familly_multi <- function() {

  markers <- c("GH","ITS","18S")

  # ---- Botanique ----
  tab_bota <- as.data.frame(table(bota$Family)) %>%
    rename(Clade = Var1, Bota_ab_rel = Freq) %>%
    mutate(Bota_ab_rel = Bota_ab_rel * 100 / sum(Bota_ab_rel))

  # ---- Génétique pour tous les marqueurs ----
  df_all <- map_dfr(markers, function(marker){

    tab_genet <- tab_rel_ab(Marker = marker) %>%
      mutate(
        Genet_ab_rel = as.numeric(Genet_ab_rel),
        Marker = marker
      )

    full_join(tab_genet, tab_bota, by = "Clade")

  })

  # ---- Corrélations ----
  corr_df <- df_all %>%
    group_by(Marker) %>%
    summarise(
      r = cor(Genet_ab_rel, Bota_ab_rel, use = "complete.obs"),
      p = cor.test(Genet_ab_rel, Bota_ab_rel)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      r = round(r,3),
      p_text = case_when(
        p < 0.001 ~ "< 0.001",
        p < 0.01  ~ "< 0.01",
        p < 0.05  ~ "< 0.05",
        p < 0.1   ~ "< 0.1",
        TRUE ~ "n.s."
      )
    )

  # ---- Format large pour les fonctions percent ----
  df_wide <- df_all %>%
    select(Clade, Marker, Genet_ab_rel, Bota_ab_rel) %>%
    pivot_wider(names_from = Marker, values_from = Genet_ab_rel)

  # ---- Out of bota ----
  df_wide <- percent_out_bota(
    df = df_wide,
    var_ref = "Clade",
    var_per_cols = markers,
    var_na = "Bota_ab_rel"
  )

  # ---- Groupement <1% ----
  df_wide <- percent_other(
    df = df_wide,
    var_ref = "Clade",
    var_per_cols = c("Bota_ab_rel", markers),
    threshold = 1
  )

  # ---- Format long pour ggplot ----
  df_long <- df_wide %>%
    pivot_longer(
      cols = c(Bota_ab_rel, GH, ITS, `18S`),
      names_to = "Label",
      values_to = "Abondance"
    ) %>%
    mutate(
      Label = recode(Label,
                     "Bota_ab_rel" = "Bota"),
      Label = factor(Label, levels = c("Bota","GH","18S","ITS"))
    ) %>%
    complete(
      Label,
      Clade,
      fill = list(Abondance = 0)
    )

  # ---- Couleurs ----
  my_colors <- c(
    "Bota" = "#33a02c",
    "GH"   = "#1f78b4",
    "ITS"  = "#ff7f00",
    "18S"  = "#e31a1c"
  )

  # ---- Barplot ----
  p <- ggplot(df_long, aes(x = Abondance, y = reorder(Clade, Abondance), fill = Label)) +
    geom_col(position = "dodge") +
    labs(
      x = "Abondance relative (%)",
      y = "Clade",
      fill = "Type/Marqueur",
      title = "Abondance relative par famille selon l'inventaire Botanique et Génétique"
    ) +
    scale_fill_manual(values = my_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      text = element_text(size = 13),
      axis.line = element_line(color = "gray"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    )

  # ---- Annotation corrélations ----
  for(i in seq_len(nrow(corr_df))){

    p <- p + annotate(
      "text",
      x = max(df_long$Abondance) * 0.7,
      y = i*3 + 3,
      label = paste0(
        "Pearson ", corr_df$Marker[i],
        ":\nr = ", corr_df$r[i],
        "\np ", corr_df$p_text[i]
      ),
      hjust = 0.5,
      size = 4
    )

  }

  return(p)
}

plot_compare_OTU_grid <- function(){
  
  grid_pred <- genet_grid_read(
    PROJECT_PATH,
    "GH", # Marker
    "Family",  # Level
    "Fabaceae",    # Name
    total = FALSE
    ) %>%
    mutate(OTU = replace_na(OTU, 0))
  
  grid_obs  <- bota_grid_compute(
    bota_df    = bota,
    level      = "Family",
    name       = "Fabaceae",
    var_bota   = "OTU",
    user_norm  = FALSE,
    radius = 10
    )
  
  plot_genet <- plot_grid(grid_pred,var="OTU")+ labs(title = "Genet OTU Richness Fabaceae")
  plot_bota <- plot_grid(grid_obs,var="OTU")+ labs(title = "Bota Richness Fabaceae (R=10)")
  plot_bota | plot_genet 
  
}

plot_compare_div_grid <- function(){
  
  grid_pred <- genet_grid_read(
    PROJECT_PATH,
    "GH", # Marker
    "Family",  # Level
    "Fabaceae",    # Name
    total = FALSE
  )
  
  grid_obs  <- bota_grid_compute(
    bota_df    = bota,
    level      = "Family",
    name       = "Fabaceae",
    var_bota   = "Shannon",
    user_norm  = FALSE,
    radius = 15
  )
  
  plot_genet <- plot_grid(grid_pred,var="Shannon") + 
    labs(title = "Genet Shannon Fabaceae")
  plot_bota <- plot_grid(grid_obs,var="Shannon") + 
    labs(title = "Bota Shannon Fabaceae (R=15)")
  
  plot_bota | plot_genet 
  
}


barplot_class_rank <- function(Marker = "GH"){
  library(ggplot2)
  library(ggrepel)
  
  # Liste des marqueurs
  markers <- c("GH", "ITS", "18S")
  df_list <- list()
  
  # ---- Boucle sur chaque marqueur ----
  for(marker in markers) {
    
    genet <- select_genet(gene, Marker = marker) %>%
      group_by(id_rank) %>%
      summarise(
        rank_count = n()
      ) %>%
      mutate(
        rank_percent = rank_count*100 / sum(rank_count)
      )
    
    genet_percent <- percent_other(
      df = genet,
      var_ref = "id_rank",
      var_per = "rank_percent",
      treshold = 5
    ) %>%
      arrange(desc(rank_percent)) %>%
      mutate(id_rank = factor(id_rank, levels = id_rank))
    
    # Ajouter colonne Marker
    df_rel_ab$Marker <- marker
    df_list[[marker]] <- df_rel_ab
  }
  
  # Fusionner tous les marqueurs
  df_all <- bind_rows(df_list)
  
  # Mettre en format long
  df_long <- df_all %>%
    pivot_longer(cols = c(Genet_ab_rel, Bota_ab_rel),
                 names_to = "Type",
                 values_to = "Abondance") %>%
    mutate(
      Label = case_when(
        Type == "Bota_ab_rel" ~ "Bota",
        Marker == "GH" ~ "GH",
        Marker == "ITS" ~ "ITS",
        Marker == "18S" ~ "18S"
      ),
      Label = factor(Label, levels = c("GH", "18S", "ITS"))  # ordre vertical des barres
    )
  
  # Rank orader
  rank_order <- c("cellular root","domain","kingdom","phylum","subphylum","clade","class","order","family","subfamily","tribe","genus","subgenus","species","subspecies","varietas")
  
  # Couleurs
  my_colors <- c("GH" = "#1f78b4", "ITS" = "#ff7f00", "18S" = "#e31a1c")
  
  # Barplot
  ggplot(df_long, aes(x = Abondance, y = reorder(Clade, Abondance), fill = Label)) +
    geom_col(position = "dodge") +
    labs(
      x = "Pourcentage des calssification (%)",
      y = "Rank",
      fill = "Marqueur",
      title = "Précision de la classification en selon le marqueur génétique"
    ) +
    scale_fill_manual(values = my_colors)

  
}

barplot_class_rank <- function(){
  
  library(dplyr)
  library(ggplot2)
  library(forcats)
  
  markers <- c("GH","18S","ITS")
  df_list <- list()
  
  for(marker in markers){
    
    genet <- select_genet(gene, Marker = marker) %>%
      group_by(id_rank) %>%
      summarise(rank_count = n(), .groups = "drop") %>%
      mutate(rank_percent = rank_count*100/sum(rank_count))
    
    genet_percent <- percent_other(
      df = genet,
      var_ref = "id_rank",
      var_per = "rank_percent",
      treshold = 5
    )
    
    genet_percent$Marker <- marker
    df_list[[marker]] <- genet_percent
  }
  
  df_all <- bind_rows(df_list)
  
  rank_order <- c(
    "Other","cellular root","domain","kingdom","phylum","subphylum","clade",
    "class","order","family","subfamily","tribe","genus",
    "subgenus","species","subspecies","varietas"
  )
  
  df_all <- df_all %>%
    mutate(
      id_rank = factor(id_rank, levels = rank_order),
      Marker = factor(Marker, levels = c("GH","18S","ITS"))
    )
  
  df_all_complete <- df_all %>%
    mutate(id_rank = fct_drop(id_rank)) %>%   # enlève les niveaux absents
    complete(
      Marker,
      id_rank,
      fill = list(
        rank_count = 0,
        rank_percent = 0
      )
    )
  
  ggplot(df_all_complete, aes(x = rank_percent, y = id_rank, fill = Marker)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      colour = "black",
      linewidth = 0.2
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      text = element_text(size = 14)
    ) +
    labs(
      x="Percent of classified OTU (%)",
      y="Maximal classification rank"
    )
  
}

# percent_other <- function(df, var_ref, var_per, treshold){
#   # Sum percents < treshold
#   other_percent <- df %>%
#     filter(.data[[var_per]] < treshold) %>%
#     summarise(sum = sum(.data[[var_per]])) %>%
#     pull(sum)
#   # Create line "Other"
#   line_other <- data.frame(
#     value_ref = paste0("<",treshold,"%"),
#     value_per = other_percent
#   )
#   # Names
#   names(line_other) <- c(var_ref, var_per)
#   # Join
#   df <- df %>%
#     filter(.data[[var_per]] >= treshold) %>% # remove lines < treshold
#     bind_rows(line_other)        # Add line sum "Other"
#   
#   return(df)
# }

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

percent_out_bota <- function(df, 
                             var_ref ="Clade", 
                             var_per = "Genet_ab_rel", 
                             var_na = "Bota_ab_rel"){
  
  # Somme des pourcentages hors bota
  other_percent <- df %>%
    dplyr::filter(is.na(.data[[var_na]])) %>%
    dplyr::summarise(sum = sum(.data[[var_per]], na.rm = TRUE)) %>%
    dplyr::pull(sum)
  
  # Ligne Out_of_bota
  line_other <- tibble::tibble(
    !!var_ref := "Out_of_bota",
    !!var_per := other_percent
  )
  
  # Conserver uniquement les lignes avec Bota
  df <- df %>%
    dplyr::filter(!is.na(.data[[var_na]])) %>%
    dplyr::bind_rows(line_other)
  
  return(df)
}

multi_filter_flow <- function(){
  
  # Init
  markers <- c("GH","18S","ITS")
  df_list <- list()
  # Loop on markers
  for(marker in markers){
    
    # marker taxid
    taxid_vec <- unique(gene[[marker]]$taxid)
    
    # filter tax_table
    df <- tax_table %>%
      filter(taxid %in% taxid_vec)
    
    # Filtering steps dataframe
    filter_steps <- tibble(
      step = c(
        "All species",
        "Plants only",
        "Dicots",
        "Family identified",
        "Species identified",
        "Valid sp name",
        "In Bota Inventory"
      ),
      Marker = marker,
      n = c(
        nrow(df),
        sum(df$is_plant, na.rm = TRUE),
        sum(df$is_plant & !df$is_monocot, na.rm = TRUE),
        sum(df$is_plant & !df$is_monocot & df$not_na_family, na.rm = TRUE),
        sum(df$is_plant & 
              !df$is_monocot &
              df$not_na_family & 
              df$not_na_sp, 
            na.rm = TRUE),
        sum(df$is_plant & 
              !df$is_monocot &
              df$not_na_family & 
              df$not_na_sp & 
              df$valid_sp_name, 
            na.rm = TRUE),
        sum(df$is_plant & 
              df$not_na_family & 
              df$not_na_sp & 
              df$valid_sp_name &
              df$is_in_bota, 
            na.rm = TRUE)
      )
    ) %>% 
      mutate(`%` = (n*100) / nrow(df)) %>%
      mutate(step = fct_reorder(step, n, .desc = TRUE)) %>%
      mutate(Marker = factor(Marker, levels = c("GH", "18S", "ITS")))
    # Store in df_list
    df_list[[marker]] <- filter_steps
  }
  
  df_all <- bind_rows(df_list)
  
  ggplot(df_all, aes(x = step, y = n, fill = Marker)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      colour = "black",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = n), 
      vjust = -0.5,
      position = position_dodge(width = 0.8)
      ) +
    theme_minimal() +
    labs(
      x = "Filtering step",
      y = "Number of OTU",
      title = "OTU filtering workflow"
    ) +
    theme(text = element_text(size = 14)) + 
    geom_hline(
      yintercept = length(unique(bota$Species)),
      linetype = "dashed"
    )
  
}

box_method_corr <- function(){
  
  # ---- Species ----
  # filter GH
  taxid_vec <- unique(gene[["GH"]]$taxid)
  df <- tax_table %>%
    filter(taxid %in% taxid_vec)
  
  sp_focus <- df %>% filter(df$is_plant & 
                              df$not_na_family & 
                              df$not_na_sp & 
                              df$valid_sp_name &
                              df$is_in_bota, 
                            na.rm = TRUE)
  
  sp_vec <- sp_focus$species
  
  # Pour chaque espèce (sp_name)
  
  # ---- Genet ----
  gene_df <- genet_grid_read(m = "GH",
                             l = "Species",
                             n = "Croton matourensis", 
                             total = FALSE) %>% 
    select(x, y, Weight)
  
  # ---- Bota ----
  radii <- seq(5, 100, 5)
  
  bota_rad <- map(radii, function(R){
    
    bota_grid_compute(
      bota_df   = bota,
      name      = "Croton matourensis",
      var_bota  = "Count",
      user_norm = FALSE,
      radius    = R
    ) %>%
      mutate(method = as.character(R))
    
  })
  
  bota_rad <- bind_rows(bota_rad)
  
  # ---- Join ----
  # Ajouter DNA
  gene_long <- gene_df %>%
    rename(value = NCI) %>%
    mutate(method = "NCI")
  
  # Fusion long
  final_long <- bind_rows(gene_long, bota_rad)
  
  # Garcer l'ordre
  final_long$method <- factor(
    final_long$method,
    levels = c("DNA", as.character(radii))
  )
  
  # ---- Boxplot ----
  ggplot(final_long, aes(x = method, y = value)) +
    geom_boxplot() +
    labs(x = "Radius / Method", y = "OTU richness") +
    theme_bw()
  
  
}

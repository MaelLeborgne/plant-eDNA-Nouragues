
# ===== plots.R =====
# Regroup plot functions by Subject

# ===== Libraries =====
library(tidyverse) # Basic
library(ggplot2) # Plot
library(viridis) # turbo
library(patchwork)
library(akima)  # pour l'interpolation 2D

# ===== < UPDATED > =====

# ==== BOTA ====

# Return the map of stems
plot_bota <- function(Level){
  
  # Configuration
  config <- list(Year=2012)
  bota <- import_bota(config)
  
  # Create plot
  p <- ggplot(bota, aes(x = Xfield, y = Yfield)) +
    geom_point(
      aes(color = .data[[Level]], size = CircCorr),
      shape = 21,       # Cercle (21–25 permettent un remplissage)
      fill = NA,        # Creux (pas de remplissage)
      stroke = 1        # Épaisseur du contour
    ) +
    theme_minimal() +
    coord_fixed() +
    theme(legend.position = "none")
  
  return(p)
}

plot_grid <- function(grid, stems = FALSE, bota_focus = NULL) {
  
  # Configuration
  df <- grid$df
  config <- grid$config
  
  # Création du ggplot de base
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_tile(aes(fill = Metric)) +
    coord_fixed() +
    theme_minimal() +
    labs(fill = config$Metric ) +
    scale_fill_viridis_c(option = "turbo")
  
  # Ajouter les stems si demandé
  if(stems) {
    if(is.null(bota_focus)) stop("Vous devez fournir bota_f si stems = TRUE")
    p <- p +
      geom_point(
        data = bota_focus,
        aes(x = Xfield, y = Yfield, size = CircCorr/10),
        shape = 21,   # cercle avec contour
        fill = NA,    # pas de remplissage
        colour = "white",
        stroke = 0.3
      )
  }
  
  return(p)
}

# Fonction d'interpolation 2D
interpolate_grid_2D <- function(grid = genet_grid_read(), factor = 3) {
  
  # Configuration
  G <- grid
  df <- G$df
  
  # 1 Grille fine
  dx <- min(diff(sort(unique(df$x)))) / factor
  dy <- min(diff(sort(unique(df$y)))) / factor
  x_fine <- seq(min(df$x), max(df$x), by = dx)
  y_fine <- seq(min(df$y), max(df$y), by = dy)
  
  # 2 Interpolation 2D
  interp_res <- akima::interp(
    x = df$x,
    y = df$y,
    z = df$Metric,
    xo = x_fine,
    yo = y_fine,
    linear = FALSE, # spline bicubique
    extrap = TRUE
  )
  
  # 3 Transformer en data.frame pour ggplot
  interp_df <- expand.grid(x = interp_res$x, y = interp_res$y)
  interp_df$krig <- as.vector(interp_res$z)
  
  # 4 Tronquer à 0
  interp_df$krig <- pmax(as.vector(interp_res$z), 0)
  
  return(interp_df)
}

# Fonction de plot
plot_grid_interp <- function(grid, factor = 2) {
  
  fine <- interpolate_grid_2D(grid, factor)
  
  # Valeurs limites pour positionner l'échelle
  x_min <- min(fine$x)
  x_max <- max(fine$x)
  y_min <- min(fine$y)
  y_max <- max(fine$y)
  
  n <- 256
  
  ggplot(fine, aes(x = x, y = y, fill = krig)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = turbo(n, alpha = 1, begin = 0, end = 1, direction = 1)
    ) +
    coord_fixed() +
    theme_minimal()+
    
    # supprimer axes et ticks
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",  # légende sous le plot
      legend.margin = margin(t = -10),
      plot.caption = element_text(
        hjust = 0.5,   # centré
        size = 14
      )
    ) +
    
    # légende horizontale sans titre
    guides(
      fill = guide_colorbar(
        title = NULL, 
        direction = "horizontal",
        barwidth = unit(8, "lines"),  # "npc" = 1 fois la largeur du panel
        barheight = unit(0.5, "lines")  # hauteur de la barre
      )
    )+
    
    # ligne d'échelle horizontale (50 unités)
    geom_segment(
      aes(
        x = x_min + 20,
        xend = x_min + 20 + 50,
        y = y_max - 0.05 * (y_max - y_min),
        yend = y_max - 0.05 * (y_max - y_min)
      ),
      color = "white",
      size = 0.7
    ) +
    
    # label "50 m"
    annotate(
      "text",
      x = x_min + 20 + 25,  # au milieu de la ligne
      y = y_max - 0.08 * (y_max - y_min),
      label = "50 m",
      hjust = 0.5,
      vjust = 1,
      size = 4,
      color = "white"
    )+
    # ajouter nom de la famille centré en bas
    labs(caption = grid$config$Name)
}

# ==== PATCH ====

patch_turbo_family <- function(factor = 1){
  
  # Configuration
  basic_config <- list(
    Type = "genet",
    Rar = FALSE,
    Marker = "GH",
    Metric = "OTU",
    Level = "Family"
  )
  
  # indices lettres pour sous-plots
  letters_idx <- letters[1:26]  # a-z
  
  # 13 familles les plus abondante (>1%, voir graph plot_tab_familly_multi())
  familles <- c("Fabaceae","Sapotaceae","Lecythidaceae","Moraceae","Chrysobalanaceae","Annonaceae","Putranjivaceae","Meliaceae","Apocynaceae","Malvaceae","Urticaceae","Euphorbiaceae","Sapindaceae")
  
  plots <- list()
  
  for(i in seq_along(familles)) {
    
    # Configuration
    fam <- familles[i]
    config <- basic_config
    config$Name <- familles[i]
    
    # Test
    print(config$Name)
    
    # Read grid
    grid <- grid_read(config)
    
    # générer la heatmap
    p <- plot_grid_interp(grid, factor = factor)
    
    # ajouter indice en haut à gauche
    p <- p +
      annotate(
        "text",
        x = min(grid$x),
        y = max(grid$y + 50),
        label = paste0("(", letters_idx[i], ")"),
        fontface = "bold",
        hjust = 0,
        vjust = 1,
        size = 6,
        color = "black"
      )
    
    plots[[i]] <- p
  }
  
  # assembler les plots en patchwork
  patchwork::wrap_plots(plots, ncol = 7)
}

ab_rel_plot <- function(df_plot){
  
  # s'assurer que Clade est un facteur dans l'ordre de remplissage désiré
  # par exemple, ici on garde l'ordre alphabétique pour le fill
  df_plot_rar <- df_plot_rar %>%
    mutate(Clade = factor(Clade, levels = sort(unique(Clade))))
  
  # calculer les positions exactes pour les labels
  df_labels <- df_plot_rar %>%
    group_by(Source) %>%
    # ordre du facteur Clade pour correspondre à l'empilement ggplot2
    arrange(Source, desc(Clade)) %>%
    mutate(ymax = cumsum(Metric),
           ymin = ymax - Metric,
           ymid = (ymin + ymax)/2) %>%
    ungroup() %>%
    filter(Metric > 2) %>%  # labels seulement pour segments >5%
    mutate(Label = Clade)
  
  # plot
  ggplot(df_plot_rar, aes(x = Source, y = Metric, fill = Clade)) + 
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
    
    geom_label(
      data = df_labels,
      aes(y = ymid, label = Label),
      fill = "black",         # fond noir
      color = "white",        # texte blanc
      fontface = "bold.italic",
      label.size = 0,
      size = 3
    ) +
    
    scale_fill_viridis(discrete = TRUE, option = "C") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      x = "Source",
      y = "Proportion (%)",
      fill = "Clade",
      title = "Proportion par Clade",
      subtitle = "Abondances relatives par Source"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}



# ===== 2. CORRELATION TABLE =====


corr_tab <- function(config_tab){
  
  # ordre des colonnes (Metric)
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("Count","NCI","T+","T-","MCC")
  )
  
  # ordre des lignes (Level)
  config_tab$Level <- factor(
    config_tab$Level,
    levels = c("Family","Genus","Species")
  )
  
  # ordre alphabétique des Name
  config_tab$Name <- factor(
    config_tab$Name,
    levels = rev(sort(unique(config_tab$Name)))
  )
  
  ggplot(config_tab,
         aes(x = Metric_corr,
             y = Name,
             fill = Correlation)) +
    
    geom_tile(color = "grey90", linewidth = 0.3) +
    
    geom_text(aes(label = sprintf("%.2f", Correlation)),
              size = 3) +
    
    facet_grid(Level ~ g.Marker,
               scales = "free",
               space = "free") +
    
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "white"
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold")
    )
}

# ===== 3. CORRELOGRAME RADIUS =====

plot_correlogramme <- function(tab_rad, Marker){
  tab_summary <- tab_rad %>%
    filter(g.Marker == Marker) %>%
    group_by(b.Radius, Level) %>%
    summarise(mean_cor = mean(Correlation, na.rm = TRUE),
              sd_cor = sd(Correlation, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(tab_summary, aes(x = b.Radius, y = mean_cor, color = Level, group = Level)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_cor - sd_cor, ymax = mean_cor + sd_cor), width = 0.1) +
    labs(x = "Radius (m)", y = "Correlation") +
    theme_bw() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) + 
    theme(plot.title = element_text(
      hjust = 0.5,   # centré
      size = 14,
      face = "bold"
    )) +
    labs(title = Marker)
}




# ===== < DEPRECATED > =====

mono_turbo_species <- function(grid, bota, factor = 2) {
  
  fine <- plot_grid_interp(grid, factor)
  
  # Valeurs limites pour positionner l'échelle
  x_min <- min(fine$x)
  x_max <- max(fine$x)
  y_min <- min(fine$y)
  y_max <- max(fine$y)
  
  n <- 256
  
  ggplot(fine, aes(x = x, y = y, fill = krig)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = turbo(n, alpha = 1, begin = 0, end = 1, direction = 1),
      breaks = scales::breaks_pretty(n = 2)
    ) +
    coord_fixed() +
    theme_minimal()+
    
    # supprimer axes et ticks
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",  # légende sous le plot
      legend.margin = margin(t = -10),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,   # centré
        size = 9,
        face = "bold.italic"
      )
    ) +
    
    # légende horizontale sans titre
    guides(
      fill = guide_colorbar(
        title = NULL, 
        direction = "horizontal",
        barwidth = unit(5, "lines"),  # "npc" = 1 fois la largeur du panel
        barheight = unit(0.5, "lines")  # hauteur de la barre
      )
    ) +
    
    # ligne d'échelle horizontale (50 unités)
    geom_segment(
      aes(
        x = x_min + 20,
        xend = x_min + 20 + 50,
        y = y_max - 0.05 * (y_max - y_min),
        yend = y_max - 0.05 * (y_max - y_min)
      ),
      color = "white",
      size = 0.5
    ) +
    
    # label "50 m"
    annotate(
      "text",
      x = x_min + 20 + 25,  # au milieu de la ligne
      y = y_max - 0.08 * (y_max - y_min),
      label = "50 m",
      hjust = 0.5,
      vjust = 1,
      size = 3,
      color = "white"
    )  +
    geom_point(
      data = bota,
      aes(x = Xfield, y = Yfield, size = CircCorr/10),
      shape = 21,   # cercle avec contour
      fill = NA,    # pas de remplissage
      colour = "white",
      stroke = 0.3,
      show.legend = FALSE
    ) +
    # ajouter nom de la famille centré en bas
    labs(caption = grid$config$Name)
}


patch_turbo_species <- function(){
  
  genet_sp <- give_genet_taxa(Rar = T, Level = "Species", Marker = "GH")
  bota_sp <- import_bota(list("Year"=2022)) %>%
    pull(Species) %>%
    unique()
  common_sp <- intersect(genet_sp, bota_sp)
  
  plots <- list()
  
  for(i in seq_along(common_sp)) {
    sp <- common_sp[i]
    
    config <- list(
      "Type" = "genet",
      "Rar" = TRUE,
      "Marker" = "GH",
      "Level" = "Species",
      "Name" = sp,
      "Metric" = "Reads"
    )
    
    genet_grid <- grid_read(config)
    
    bota_focus <- import_bota(list("Year"=2022)) %>%
      filter(Species == sp)
    
    print(sp)
    
    # générer la heatmap
    p <- mono_turbo_species(genet_grid, bota_focus)
    
    # ajouter indice en haut à gauche
    p <- p +
      annotate(
        "text",
        x = min(genet_grid$x),
        y = max(genet_grid$y + 60),
        label = paste0("(", i, ")"),
        fontface = "bold",
        hjust = 0,
        vjust = 1,
        size = 5,
        color = "black"
      )
    
    plots[[i]] <- p
  }
  
  # assembler les plots en patchwork
  patchwork::wrap_plots(plots, ncol = 12)
  
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

# ==== PATCH ====

plot_corr_heatmap_families <- function(){
  
  # mêmes familles que ton patchwork
  familles <- c(
    "Fabaceae","Sapotaceae","Lecythidaceae","Moraceae",
    "Chrysobalanaceae","Annonaceae","Putranjivaceae","Meliaceae",
    "Apocynaceae","Malvaceae","Urticaceae","Euphorbiaceae","Sapindaceae"
  )
  
  # récupérer les grilles
  grids <- lapply(familles, function(fam){
    
    g <- genet_grid_read(
      m = "GH",
      l = "Family",
      n = fam,
      from_csv = TRUE,
      total = FALSE
    )
    
    g %>%
      select(x,y,OTU) %>%
      mutate(Family = fam)
    
  })
  
  grid_all <- bind_rows(grids)
  
  # transformer en matrice
  mat <- grid_all %>%
    unite(cell, x, y, remove = FALSE) %>%
    select(cell, Family, OTU) %>%
    pivot_wider(names_from = Family, values_from = OTU)
  
  mat_values <- mat %>% select(-cell)
  
  # corrélation
  corr_mat <- cor(mat_values, use = "pairwise.complete.obs", method = "pearson")
  
  corr_df <- as.data.frame(as.table(corr_mat))
  colnames(corr_df) <- c("Family1","Family2","Correlation")
  
  # heatmap
  p <- ggplot(corr_df, aes(Family1, Family2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-1,1),
      breaks = c(-1,-0.5,0,0.5,1)
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(
      title = "Corrélation spatiale entre familles (heatmaps GH)",
      fill = "Pearson r"
    )
  
  return(p)
}





patch_compare <- function(var_bota,
                          var_genet,
                          b_grid = bota_grid(),
                          g_grid = genet_grid(),
                          b_f    = bota_f(),
                          stems  = FALSE,
                          Title = "") {
  
  # --- Plot 1 : Bota density ---
  p1 <- plot_grid(
    grid  = b_grid,
    var   = var_bota,
    stems = stems,
    bota_f = b_f
  ) +
    labs(title = "Bota (R=10)")
  
  # --- Plot 2 : Genet density ---
  p2 <- plot_grid(
    grid  = g_grid,
    var   = var_genet,
    stems = stems,
    bota_f = b_f
  ) +
    labs(title = "Genet")
  
  df_compare <- df_compare_join(grid_pred  = g_grid, 
                                grid_obs   = b_grid, 
                                var_pred   = var_genet, 
                                var_obs    = var_bota,
                                total = FALSE)
  
  p3 <- plot_compare(df_compare, bota_f = b_f, stems = stems)
  
  # Combine en ligne avec patchwork
  combined_plot <- ( p1 | p2 | p3 ) +
    plot_annotation(
      title = Title
    )
  
  return(combined_plot)
}

# Return the histogram of the grid (bota or genet)
plot_hist <- function(grid, col_name){
  # Create the ggplot histogram
  p <- ggplot(grid, aes(x = .data[[col_name]])) +
    geom_histogram() +
    theme_minimal()
  # Return
  return(p)
}



# ===== 1. FAMILY COMPOSITION BY MARKER =====






# ggplot(df, aes(x = group, fill = category)) +
#   geom_bar(position = "fill")

plot_tab_familly_multi <- function() {
  
  markers <- c("GH","ITS","18S")
  
  # Botanique
  tab_bota <- as.data.frame(table(bota$Family)) %>%
    rename(Clade = Var1, Bota_ab_rel = Freq) %>%
    mutate(Bota_ab_rel = Bota_ab_rel * 100 / sum(Bota_ab_rel))
  
  # Génétique pour tous les marqueurs
  df_all <- map_dfr(markers, function(marker){
    
    tab_genet <- tab_rel_ab(Marker = marker) %>%
      mutate(
        Genet_ab_rel = as.numeric(Genet_ab_rel),
        Marker = marker
      )
    
    full_join(tab_genet, tab_bota, by = "Clade")
    
  })
  
  # Corrélations
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
  
  # Format large pour les fonctions percent
  df_wide <- df_all %>%
    select(Clade, Marker, Genet_ab_rel, Bota_ab_rel) %>%
    pivot_wider(names_from = Marker, values_from = Genet_ab_rel)
  
  # Out of bota
  df_wide <- percent_out_bota(
    df = df_wide,
    var_ref = "Clade",
    var_per_cols = markers,
    var_na = "Bota_ab_rel"
  )
  
  # Groupement <1%
  df_wide <- percent_other(
    df = df_wide,
    var_ref = "Clade",
    var_per_cols = c("Bota_ab_rel", markers),
    threshold = 1
  )
  
  # Format long pour ggplot
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
  
  # Couleurs
  my_colors <- c(
    "Bota" = "#33a02c",
    "GH"   = "#1f78b4",
    "ITS"  = "#ff7f00",
    "18S"  = "#e31a1c"
  )
  
  # Barplot
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
  
  # Annotation corrélations
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








# ==== GENET ====

# ==== COMPARE ====
# Return plot Compare
plot_compare <- function(df_compare, bota_f, stems = FALSE) {
  
  # Base plot
  p <- ggplot(df_compare, aes(x = x, y = y)) +
    geom_tile(aes(fill = err_abs)) +
    scale_fill_gradient2(
      low = "red",      # valeurs négatives
      mid = "white",    # neutre (0)
      high = "green",   # valeurs positives
      midpoint = 0
    ) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "err = Genet - Bota")
  
  # Ajouter les stems si demandé
  if(stems) {
    if(is.null(bota_f)) stop("bota_f doit être fourni si stems = TRUE")
    p <- p +
      geom_point(
        data = bota_f,
        aes(x = Xfield, y = Yfield, size = CircCorr/10),
        shape = 21,  # cercle avec contour
        fill = NA,   # creux
        colour = "black",
        stroke = 0.3
      )
  }
  
  return(p)
}

# Plot 2D Genet / parameter
plot_x_y <- function(df_compare, 
                     x_param = "Count_dens", 
                     y_param = "Weight_dens"
                     ){
  ggplot(df_compare,aes(x = value_obs, y = value_pred)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "grey70")+
    theme_minimal() +
    theme(
      axis.line = element_line(color = "gray"),  # Lignes des axes en noir
      axis.text = element_text(color = "black"),  # Texte des axes en noir
      axis.ticks = element_line(color = "black")  # Ticks des axes en noir
    ) +
    labs(
      x = x_param,
      y = y_param
    )
}


# Barplot of engin distance traveled
barplot_R2_taxon <- function(df){
  # Trier les données
  df <- df |> 
    arrange(desc(R2)) |>
    mutate(Name = factor(Name, levels = Name))
  # Plot
  ggplot(df, aes(reorder(Name, R2), R2)) +
    geom_col() +
    geom_text( # Etiquettes heures
      aes(label = paste0("R2 = ",round(R2))),
      hjust = -0.2,          # décale vers la droite
      size = 4,
      fontface="bold"
    ) +
    coord_flip() +
    expand_limits(y = max(df$R2) * 1.1) + # Laisse de la place pour les labels
    labs( # Labels
      title="Variance expliqué en fonction du taxon", 
      x="", y="R2"
    )
}

patch_BA <- function(){
  p1 <- plot_grid(bota_grid_count(),var="Count") +
    labs(title = "Bota Count (R=10)")
  p2 <- plot_grid(genet_grid_read(),var="Weight") +
    labs(title = "Genet Number of Reads")
  p4 <- plot_grid(grid_basal_area(),var="Basal_area_m2") +
    labs(title = "Bota Basal Area (R=10)")
  p1 | p2 | p4
}

plot_BA_reads <- function(){
  
  df <- df_compare_join(genet_grid_read(),grid_basal_area(),"OTU","Basal_area_m2")
  plot_x_y(df_compare = df,x_param = "value_obs", y_param = "value_pred")
}

# ---- 16 MARS MEETING ----

plot_abrel_comp <- function(){
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
  
  # Genet
  gene_df <- genet_grid_read(m = "GH", 
                             n = "All", 
                             total = FALSE) %>% 
    select(x, y, OTU) %>%
    rename(DNA = OTU)
  
  # Bota
  radii <- seq(5, 100, 5)
  
  bota_rad <- map(radii, function(R){
    
    bota_grid_read(
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
  
  # Join
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
  
  # Boxplot
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
    
    # Stem census (bota)
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
    
    # eDNA
    dna_comm <- readRDS("/home/jean/Bureau/GH_clean/GH_mblist_annot_metabar_clean_clust97_RAR.rds")
    dna_comm <- dna_comm$reads
    dna_rich <- rowSums(dna_comm > 0)
    
    # Intersect
    common <- intersect(names(stem_rich), names(dna_rich))
    
    # Pearson correlation
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
  
  # Bind
  cor_results <- bind_rows(cor_results)
  
  # Plot
  
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




plot_compare_OTU_grid <- function(){
  
  grid_pred <- genet_grid_read(
    PROJECT_PATH,
    "GH", # Marker
    "Family",  # Level
    "Fabaceae",    # Name
    total = FALSE
  ) %>%
    mutate(OTU = replace_na(OTU, 0))
  
  grid_obs  <- bota_grid_read(
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
  
  grid_obs  <- bota_grid_read(
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
  
  # Boucle sur chaque marqueur
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



box_method_corr <- function(){
  
  # Species
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
  
  # Genet
  gene_df <- genet_grid_read(m = "GH",
                             l = "Species",
                             n = "Croton matourensis", 
                             total = FALSE) %>% 
    select(x, y, Weight)
  
  # Bota
  radii <- seq(5, 100, 5)
  
  bota_rad <- map(radii, function(R){
    
    bota_grid_read(
      bota_df   = bota,
      name      = "Croton matourensis",
      var_bota  = "Count",
      user_norm = FALSE,
      radius    = R
    ) %>%
      mutate(method = as.character(R))
    
  })
  
  bota_rad <- bind_rows(bota_rad)
  
  # Join
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
  
  # Boxplot
  ggplot(final_long, aes(x = method, y = value)) +
    geom_boxplot() +
    labs(x = "Radius / Method", y = "OTU richness") +
    theme_bw()
  
}

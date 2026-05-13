
# ===== plots.R =====
# Regroup plot functions by Subject

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # select, filter
library(ggplot2) # Plot
library(viridis) # turbo
library(patchwork)
library(akima)  # pour l'interpolation 2D
library(emmeans)
library(multcomp)
library(ggtext)

# ==== 0. BOTA ====

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



# ==== 1. RELATIVE ABUNDANCE ====

ab_rel_tab <- function(df_ab_rel){
  
  # Ordre des Sources
  df_ab_rel$Source <- factor(
    df_ab_rel$Source,
    levels = c("Bota","GH","ITS","18S")
  )
  
  # Ordre des Clades
  order_bota <- df_ab_rel %>%
    filter(Source == "Bota") %>%
    arrange(desc(Metric)) %>%
    distinct(Clade) %>%
    pull(Clade)
  
  # autres Clade (non présents dans Bota)
  other_clade <- setdiff(unique(df_ab_rel$Clade), order_bota)
  
  # combiner les deux
  order_clade <- c(order_bota, other_clade)
  
  # Appliquer l'ordre
  df_ab_rel$Clade <- factor(df_ab_rel$Clade, levels = rev(order_clade))
  
  # Clean data
  tab_clean <- droplevels(na.omit(df_ab_rel))
  
  # Labels conditionels Y
  labels_y <- tab_clean %>%
    group_by(Clade) %>%
    summarise(Metric = max(Metric), .groups = "drop") %>%
    mutate(
      label = if_else(
        Metric > 2,
        paste0("<b>", Clade, "</b>"),
        as.character(Clade)
      )
    )
  
  # Factors
  lab_vec <- setNames(labels_y$label, labels_y$Clade)
  tab_clean$Clade <- factor(tab_clean$Clade, levels = names(lab_vec))
  
  # Plot
  ggplot(tab_clean,
         aes(x = Source,
             y = Clade,
             fill = Metric)) +
    
    geom_tile(color = "grey90", linewidth = 0.3) +
    
    geom_text(
      aes(
        label = sprintf("%.1f", Metric),
        fontface = ifelse(Metric > 2, "bold", "plain")
      ),
      size = 3
    ) +
    
    scale_x_discrete(position = "top") +
    scale_y_discrete(
      labels = lab_vec,
      guide = guide_axis(n.dodge = 1)
    ) +
    
    scale_fill_gradient2(
      low = "white",
      high = "white",
      midpoint = 0,
      limits = c(0, 100),
      na.value = "white"
    ) +
    
    theme_minimal(base_size = 10) +
    
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(face = "bold", size = 14),
      axis.text.y.left = ggtext::element_markdown(
        hjust = 1, color = "black")
    )
}

ab_rel_plot <- function(df_plot){
  
  # Mettre les clades dans l'ordre alphabétique
  df_plot <- df_plot %>%
    mutate(Clade = factor(Clade, levels = sort(unique(Clade))))
  
  # Choisir l'ordre des Source
  df_plot <- df_plot %>%
    mutate(Source = factor(Source, levels = c("Bota","GH","ITS","18S")))
  
  # calculer les positions exactes pour les labels
  df_labels <- df_plot %>%
    group_by(Source) %>%
    arrange(Source, desc(Clade)) %>%
    mutate(
      ymax = cumsum(Metric),
      ymin = ymax - Metric,
      ymid = (ymin + ymax)/2
    ) %>%
    ungroup() %>%
    mutate(
      Label = case_when(
        Clade %in% c("<2%","Indet.") ~ paste0("<span style='color:black;'>", Clade, "</span>"),
        in_bota                      ~ paste0("<span style='color:forestgreen;'>", Clade, "</span>"),
        in_both                      ~ paste0("<span style='color:blue;'>", Clade, "</span>"),
        in_genet                     ~ paste0("<span style='color:red;'>", Clade, "</span>"))
    )
  
  # Labels vector
  labels_vec <- setNames(df_labels$Label, df_labels$Clade)
  
  # x labels colors
  labels_source <- c(
    "Bota" = "<span style='color:forestgreen;'> Bota </span>",
    "GH" = "<span style='color:red;'> GH </span>",
    "ITS" = "<span style='color:red;'> ITS </span>",
    "18S" = "<span style='color:red;'> 18S </span>"
  )
  
  # Plot
  ggplot(df_plot, aes(x = Source, y = Metric, fill = Clade)) + 
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
    
    ggtext::geom_richtext(
      data = df_labels,
      aes(y = ymid, label = Label),
      fill = "white",
      label.size = 0,
      fontface = "bold.italic",
      size = 3
    ) +
    
    scale_fill_viridis(discrete = TRUE, option = "C") +
    scale_x_discrete(labels = labels_source) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      x = "Inventaire / Marqueur",
      y = "Proportion (%)",
      fill = "Clade",
      title = "Abondances relatives par Famille",
      subtitle = "En richesse botanique et richesse d'OTU"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
      axis.text.x.bottom = ggtext::element_markdown(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none"
    )+ 
    geom_hline(yintercept = c(25,50,75), linetype = "dashed", color="white")
}

growth_form_plot <- function(growth_form){
  
  # Linéariser le dataframe
  growth_form_long <- growth_form %>%
    pivot_longer(
      cols = -c(family, Clade, in_bota, in_both, in_genet),
      names_to = "growth_form",
      values_to = "Percent"
    ) %>%
    mutate(
      color_order = case_when(
        in_bota ~ 1,
        in_both ~ 2,
        in_genet ~ 3
      )
    )
  
  # ordre des family selon in_bota
  family_order <- growth_form_long %>%
    distinct(family, color_order) %>%
    arrange(desc(color_order), desc(family)) %>%   # color_order, puis ordre alphabétique
    pull(family)
  
  # et alphabétique
  growth_form_long$family <- factor(
    growth_form_long$family,
    levels = family_order
  )
  
  # créer labels avec couleur conditionelles
  label_df <- growth_form_long %>%
    distinct(family, color_order) %>%
    mutate(family = as.character(family)) %>%
    mutate(
      label = case_when(
        family %in% c("<2%","Indet.") ~ paste0("<span style='color:black;'>", family, "</span>"),
        color_order == 1 ~ paste0("<span style='color:forestgreen;'>", family, "</span>"),
        color_order == 2 ~ paste0("<span style='color:blue;'>", family, "</span>"),
        color_order == 3 ~ paste0("<span style='color:red;'>", family, "</span>"))
    )
  # To Named vector
  labels_vec <- setNames(label_df$label, label_df$family)
  
  # Plot
  ggplot(growth_form_long, aes(x = Percent, y = family, fill = growth_form)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    scale_y_discrete(
      labels = labels_vec,
      guide = guide_axis(n.dodge = 1)
    ) +
    labs(
      x = "Proportion (%)",
      y = "",
      fill = "Growth Form",
      title = "Forme de croissance végétale selon la famille",
      subtitle = "En richesse spécifique sur l'ensemble des nouragues"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y.left = ggtext::element_markdown(
        hjust = 1,
        face = "bold.italic"
      )
    ) + 
    geom_vline(xintercept = c(25,50,75), linetype = "dashed")
}

# ==== 2. PATCH ====

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

# Plot interpolate grid
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
      legend.margin = ggplot2::margin(t = -10),
      plot.caption = element_text(
        hjust = 0.5,   # centré
        size = 12
      ),
      plot.margin = ggplot2::margin(0, 20, 0, 20)
    ) +
    
    # légende horizontale sans titre
    guides(
      fill = guide_colorbar(
        title = NULL, 
        direction = "horizontal",
        barwidth = unit(6, "lines"),  # "npc" = 1 fois la largeur du panel
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
      size = 3,
      color = "white"
    )+
    # ajouter nom de la famille centré en bas
    labs(caption = grid$config$Name)
}

# Fonction de plot  with soft theme
plot_grid_interp_soft <- function(grid, factor = 3) {
  
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
      colours = turbo(100)
      #limits = c(0, 30) # Homogene scaling color
    ) +
    coord_fixed() +
    theme_light()+
    
    # supprimer axes et ticks
    theme(
      legend.margin = ggplot2::margin(t = -10),
      plot.title = element_text(hjust = 0,size = 14, face="bold"),
      plot.margin = ggplot2::margin(0, 20, 0, 20)
    )
}

patch_turbo_family <- function(factor = 1){
  
  # Configuration
  basic_config <- list(
    Type = "genet",
    Rar = FALSE,
    Marker = "GH",
    Metric = "OTU",
    Level = "Family"
  )
  
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
        x = min(grid$df$x),
        y = max(grid$df$y + 60),
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
  patchwork::wrap_plots(plots, ncol = 4)
}

mono_turbo_species <- function(grid, factor = 2) {
  
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
      ),
      plot.margin = margin(0, 20, 0, 20)
    ) +
    
    # légende horizontale sans titre
    guides(
      fill = guide_colorbar(
        title = NULL, 
        direction = "horizontal",
        barwidth = unit(4, "lines"),  # "npc" = 1 fois la largeur du panel
        barheight = unit(0.5, "lines")  # hauteur de la barre
      )
    ) +
    
    # ligne d'échelle horizontale (50 unités)
    geom_segment(
      aes(
        x = x_min ,
        xend = x_min + 50,
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
    ) +
    # ajouter nom de la famille centré en bas
    labs(caption = grid$config$Name)
}


patch_turbo_species <- function(){
  load("~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.df_corr_tot.RData",
       envir = .GlobalEnv)
  
  common_sp <- df_corr_tot %>%
    dplyr::filter(
      Level %in% c("Species"),
      g.Marker %in% c("GH"),
      Metric_corr %in% c("NCI")
    ) %>%
    dplyr::select(Name,Correlation) %>%
    arrange(desc(Correlation))
  
  # Select species
  common_sp <- common_sp[c(1,7,12,13),]
  
  common_sp$plots <- lapply(seq_along(common_sp$Name),function(i){
    
    sp <- common_sp$Name[[i]]
    corr <- common_sp$Correlation[[i]] %>% round(digits = 3)
    print(sp)
    
    # Import genet grid
    genet_grid <- grid_read(list(
      "Type" = "genet",
      "Rar" = TRUE,
      "Marker" = "GH",
      "Level" = "Species",
      "Name" = sp,
      "Metric" = "Reads"
    ))
    
    # Bota year 1 grid
    bota_focus1 <- import_bota(list("Year" = 2012)) %>%
      filter(Species == sp) %>%
      mutate(Year = factor(2012, levels = c(2012, 2022)))
    
    # Bota year 2 grid
    bota_focus2 <- import_bota(list("Year" = 2022)) %>%
      filter(Species == sp) %>%
      mutate(Year = factor(2022, levels = c(2012, 2022)))
    
    # générer la heatmap
    p <- mono_turbo_species(genet_grid)
    
    # Année 1
    p <- p +
      geom_point(
        data = bota_focus1,
        aes(x = Xfield, y = Yfield, size = CircCorr/2),
        color = "red",
        shape = 21,
        stroke = 0.3,
        fill = NA
      )
    
    # Année 2
    p <- p +
      geom_point(
        data = bota_focus2,
        aes(x = Xfield, y = Yfield, size = CircCorr),
        color = "orange",
        shape = 21,
        stroke = 0.3,
        fill = NA
      )+
      guides(size = "none")
    
    # ajouter indice en haut à gauche
    p <- p +
      annotate(
        "text",
        x = (min(genet_grid$df$x)+max(genet_grid$df$x))/2,
        y = max(genet_grid$df$y + 60),
        label = paste0("R = ", corr),
        fontface = "bold",
        hjust = 0.5,
        vjust = 1,
        size = 4,
        color = "black"
      )
  })
  
  # assembler les plots en patchwork
  wrap_plots(common_sp$plots, ncol = 5)
    # plot_annotation(
    #   title = Year,
    #   theme = theme(
    #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    #   )
    # )
}

patch_turbo_species_years <- function(Year1, Year2){
  
  # Import Genet species
  genet_sp <- give_genet_taxa(Rar = T, Level = "Species", Marker = "GH")
  
  # Import bota species
  bota1 <- import_bota(list("Year" = Year1))
  bota2 <- import_bota(list("Year" = Year2))
  
  # Keep unique species
  bota_sp1 <- unique(bota1$Species)
  bota_sp2 <- unique(bota2$Species)
  
  # Common species with genet
  common_sp <- Reduce(intersect, list(genet_sp, bota_sp1, bota_sp2)) %>%
    sort()
  
  # # Select species
  # common_sp <- common_sp[c(1,8,12,13)]
  
  # List of mini plots
  plots <- list()
  
  # Loop over species
  for(i in seq_along(common_sp)) {
    
    # Set species
    sp <- common_sp[i]
    
    # Print species
    print(sp)
    
    # Genet config
    genet_config <- list(
      "Type" = "genet",
      "Rar" = TRUE,
      "Marker" = "GH",
      "Level" = "Species",
      "Name" = sp,
      "Metric" = "Reads"
    )
    
    # Genèt Grid
    genet_grid <- grid_read(genet_config)
    
    # Bota year 1 grid
    bota_focus1 <- import_bota(list("Year" = Year1)) %>%
      filter(Species == sp) %>%
      mutate(Year = factor(Year1, levels = c(Year1, Year2)))
    
    # Bota year 2 grid
    bota_focus2 <- import_bota(list("Year" = Year2)) %>%
      filter(Species == sp) %>%
      mutate(Year = factor(Year2, levels = c(Year1, Year2)))
    
    # générer la heatmap
    p <- mono_turbo_species(genet_grid)
    
    # Année 1
    p <- p +
      geom_point(
        data = bota_focus1,
        aes(x = Xfield, y = Yfield, size = CircCorr/2),
        color = "red",
        shape = 21,
        stroke = 0.3,
        fill = NA
      )
    
    # Année 2
    p <- p +
      geom_point(
        data = bota_focus2,
        aes(x = Xfield, y = Yfield, size = CircCorr),
        color = "orange",
        shape = 21,
        stroke = 0.3,
        fill = NA
      ) +
      guides(size = "none")
    
    # # Year color
    # p <- p +
    #   scale_color_manual(
    #     values = setNames(
    #       c("red", "orange"),
    #       c(as.character(Year1), as.character(Year2))
    #     ),
    #     name = "Year"
    #   ) +
    #   guides(size = "none")
    
    # index
    p <- p +
      annotate(
        "text",
        x = min(genet_grid$df$x),
        y = max(genet_grid$df$y + 60),
        label = paste0("(", i, ")"),
        fontface = "bold",
        hjust = 0,
        vjust = 1,
        size = 4
      ) 
    
    
    plots[[i]] <- p
  }
  
  # assembler les plots en patchwork
  wrap_plots(plots, ncol = 4) 
  # +
  #   plot_annotation(
  #     title = "Comparaison des Inventaires botaniques de 2012 et 2022 avec l'emprunte ADNe",
  #     subtitle = "Au niveau specifique pour le marqueur GH",
  #     theme = theme(
  #       plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #       plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12)
  #     )
  #   ) 

}





# ===== 3. CORRELATION MODELS =====


corr_tab_plot <- function(config_tab){
  
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("T+","T-","MCC","Count","NCI")
  )
  
  config_tab$Level <- factor(
    config_tab$Level,
    levels = c("Family","Genus","Species")
  )
  
  config_tab$Name <- factor(
    config_tab$Name,
    levels = rev(sort(unique(config_tab$Name)))
  )
  
  config_tab$g.Marker <- factor(
    config_tab$g.Marker,
    levels = c("GH","ITS","18S")
  )
  
  # clean data
  config_tab_clean <- droplevels(na.omit(config_tab))
  
  # Plot
  ggplot(config_tab_clean,
         aes(x = Metric_corr,
             y = Name,
             fill = Correlation)) +
    
    geom_tile(color = "grey90", linewidth = 0.3) +
    
    geom_text(aes(label = sprintf("%.2f", Correlation)),
              size = 3) +
    
    facet_grid(Level ~ g.Marker,
               scales = "free",
               space = "free",
               switch = "y",
               drop = TRUE) +
    
    scale_x_discrete(position = "top") +
    
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

corr_boxplot <- function(config_tab){
  
  # ordre des colonnes (Metric)
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("NCI","Count","T+","T-","MCC")
  )
  
  # ordre des lignes (Level)
  config_tab$Level <- factor(
    config_tab$Level,
    levels = c("Family","Genus","Species")
  )
  
  # clean data
  config_tab_clean <- droplevels(na.omit(config_tab))
  
  # Test stat global
  res_aov3 <- aov(Correlation ~ Level * Metric_corr * g.Marker, 
                  data = config_tab_clean)
  
  # Moyennes ajusté par Level
  emm3 <- emmeans(res_aov3, ~ Level | Metric_corr * g.Marker)
  
  # Lettre de significativités
  cld3 <- cld(emm3, Letters = letters, alpha = 0.1) # p-value 0.1
  cld_df3 <- as.data.frame(cld3)
  
  # Position des étiquettes
  cld_df3 <- cld_df3 %>%
    left_join(
      config_tab %>%
        group_by(Metric_corr, g.Marker, Level) %>%
        summarise(
          pos = min(Correlation, na.rm = TRUE)
        ),
      by = c("Metric_corr", "g.Marker", "Level")
    )
  
  # Plot
  ggplot(config_tab_clean,
         aes(x = Correlation,
             y = Metric_corr)) +
    
    # boxplot
    geom_violin(
      width = 0.8,
      outlier.size = 1.5,
      fill = "grey80",
      color = "black"
    ) +
    
    # points (optionnel mais recommandé)
    geom_jitter(
      width = 0.15,
      size = 1,
      alpha = 0.5
    ) +
    
    # facettes = structure tableau
    facet_grid(g.Marker ~ Level,
               scales = "free_x",
               space = "free_x",
               switch = "y",
               drop = TRUE) +
    
    # axe en haut (tableau)
    scale_x_continuous(position = "top") +
    
    geom_text(data = cld_df3,
              aes(x = pos - 0.05,   # position sur axe x
                  y = Metric_corr,
                  label = .group),
              inherit.aes = FALSE,
              size = 5,
              col = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 1, linetype = "dashed", color = "black")+
    
    theme_minimal(base_size = 11) +
    
    theme(
      panel.grid = element_blank(),
      #axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x.bottom = element_text(size=12),
      axis.line = element_line(color = "gray"),  # Lignes des axes en noir
      axis.text = element_text(color = "black"),  # Texte des axes en noir
      axis.ticks = element_line(color = "black"),  # Ticks des axes en noir
      
      # strips
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(face = "bold",size=14),
      strip.text.y.left = element_text(angle = 0, face = "bold",size=14),
      
      # axes
      axis.text.x.top = element_text(face = "bold")
    ) +
    
    labs(
      x = "Correlation coefficient"
    ) + 
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(-0.15, 1)
    )
}

corr_boxplot_taxo <- function(config_tab){
  
  # # clean data
  # config_tab_clean <- droplevels(na.omit(config_tab))
  # 
  # # Test stat global
  # res_aov3 <- aov(Correlation ~ Level, 
  #                 data = config_tab_clean)
  # 
  # # Moyennes ajusté par Level
  # emm3 <- emmeans(res_aov3, ~ Level)
  # 
  # # Lettre de significativités
  # cld3 <- cld(emm3, Letters = letters, alpha = 0.1) # p-value 0.1
  # cld_df3 <- as.data.frame(cld3)
  # 
  # # Position des étiquettes
  # cld_df3 <- cld_df3 %>%
  #   left_join(
  #     config_tab %>%
  #       group_by(Level) %>%
  #       summarise(
  #         y = min(Correlation, na.rm = TRUE)
  #       ),
  #     by = c("Level")
  #   )
  
  # Plot
  ggplot(tab,
         aes(x = Level,
             y = Correlation)) +
    
    # boxplot
    geom_boxplot(
      width = 0.8,
      outlier.size = 1.5,
      fill = "grey80",
      color = "black"
    ) +
    
    # points (optionnel mais recommandé)
    geom_jitter(
      width = 0.15,
      size = 1,
      alpha = 0.5
    ) +
    
    # facettes = structure tableau
    # facet_grid(Level ~ g.Marker,
    #            scales = "free_x",
    #            space = "free_x",
    #            switch = "y",
    #            drop = TRUE) +
    
    # axe en haut (tableau)
    # scale_x_discrete(position = "top") +
    
    # geom_text(data = cld_df3,
    #           aes(x = Metric_corr,
    #               y = y - 0.05,   # décale au-dessus
    #               label = .group),
    #           inherit.aes = FALSE,
    #           size = 5,
    #           col = "red") +
    
    theme_minimal(base_size = 11) +
    
    theme(
      #panel.grid = element_blank(),
      #axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "gray"),  # Lignes des axes en noir
      axis.text = element_text(color = "black"),  # Texte des axes en noir
      axis.ticks = element_line(color = "black"),  # Ticks des axes en noir
      
      # strips
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      
      # axes
      axis.text.x.top = element_text(face = "bold")
    ) +
    
    labs(
      y = "Presence/Absence Correlation (MCC)"
    ) +
    
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(-0.15, 1)
    )
}

plot_correlogramme <- function(tab_rad, Marker){
  # As Factors
  tab_rad$b.Radius <- as.factor(tab_rad$b.Radius)
  tab_rad$g.Marker <- as.factor(tab_rad$g.Marker)
  tab_rad$Level <- as.factor(tab_rad$Level)
  
  # # Median and quartiles
  tab <- tab_rad %>%
    filter(g.Marker == Marker)
  
  # Plot
  ggplot(tab, aes(x = b.Radius, y = Correlation, color = Level, group=Level)) +
    geom_line() +
    geom_point() +
    labs(x = "Radius (m)", y = "Correlation") +
    theme_bw() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) + 
    theme(plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    )) +
    labs(title = Marker)
}


# ==== 4. CORRELATION MORPHO ====

morpho_corr_bar <- function(){
  
  # Import Morpho
  morpho <- import_plants_morpho() %>%
    mutate(Species = gsub("_", " ", Species))
  
  # Read Model
  model <- tab_rar %>%
    filter(
      Level == "Species",
      g.Marker == "GH",
      Metric_corr == "NCI"
    ) %>%
    dplyr::select(Name, Correlation) %>%
    rename(Species = Name, Model = Correlation)
  
  # Inner Join
  model_morpho_df <- inner_join(model, morpho, by = "Species") %>%
    dplyr::select(
      -Species, -Genus, -Family, -order, -subclass2, -subclass1
    )
  
  # Prepare data
  corr_vec <- cor(model_morpho_df, use = "pairwise.complete.obs")["Model", ]
  
  corr_df <- data.frame(
    Variable = names(corr_vec),
    Correlation = corr_vec
  ) %>%
    dplyr::filter(Variable != "Model") %>%
    arrange(Correlation)
  
  # Labels colored
  colors <- list(
    "forestgreen" = c("Chlorophyll_content",	"Thickness",	"Toughness",	"Leaf_Area",	"SLA",	"C",	"N",	"13C",	"Ca",	"P",	"K"),
    "orange" = c("Trunk_bark_thickness",	"Sapwood_WSG"),
    "brown" = c("WSG"),
    "brown" = c("Diameter",	"SRL",	"Tissue_density",	"SRTA",	"branchiness")
  )
  
  # Transformer la liste en vecteur nommé (paramètre -> couleur)
  param_colors <- setNames(rep(names(colors), lengths(colors)), unlist(colors))
  
  # Générer les labels HTML automatiquement
  labels_source <- setNames(
    paste0("<span style='color:", param_colors, ";'> ", names(param_colors), " </span>"),
    names(param_colors)
  )
  
  # Plot
  ggplot(corr_df, aes(x = reorder(Variable, Correlation), y = Correlation, fill = Correlation)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Corrélation des traits morphologique avec le modèle",
      subtitle = "au niveau specifique pour le marker GH",
      x = "",
      y = "Pearson r",
      fill = "r"
    ) +
    scale_x_discrete(labels = labels_source) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
      axis.text.y.left = ggtext::element_markdown(hjust = 1, face = "bold")
    )
  
}

morpho_corr_matrix <- function(){
  
  # Morpho
  morpho <- import_plants_morpho() %>%
    mutate(Species = gsub("_", " ", Species))
  # Correlations
  model <- tab_rar %>%
    filter(Level == "Species",
           g.Marker == "GH",
           Metric_corr == "NCI") %>%
    dplyr::select(Name, Correlation) %>%
    rename(Species = Name, NCI = Correlation)
  
  # Join
  df <- inner_join(model, morpho, by = "Species") %>%
    dplyr::select(-Species, -Genus, -Family, -order, -subclass2, -subclass1)
  
  # CORRELATION
  corr_mat <- cor(df, use = "pairwise.complete.obs")
  
  # Ordre basé sur Model
  order_vars <- names(sort(corr_mat["NCI", ], decreasing = TRUE))
  order_vars <- c("NCI", setdiff(order_vars, "NCI"))
  
  # LONG FORMAT
  corr_df <- as.data.frame(as.table(corr_mat))
  colnames(corr_df) <- c("Var1", "Var2", "Correlation")
  
  # TRIANGLE (upper only)
  corr_df <- corr_df %>%
    mutate(
      Var1 = factor(Var1, levels = order_vars),
      Var2 = factor(Var2, levels = order_vars)
    ) %>%
    dplyr::filter(as.numeric(Var1) <= as.numeric(Var2))
  
  # GROUP COLORS
  colors <- list(
    forestgreen = c("Chlorophyll_content","Thickness","Toughness","Leaf_Area","SLA","C","N","13C","Ca","P","K"),
    orange = c("Trunk_bark_thickness","Sapwood_WSG"),
    brown = c("WSG","Diameter","SRL","Tissue_density","SRTA","branchiness"),
    black = c("NCI")
  )
  
  param_colors <- setNames(rep(names(colors), lengths(colors)), unlist(colors))
  
  labels <- setNames(
    paste0("<span style='color:", param_colors, ";'>", names(param_colors), "</span>"),
    names(param_colors)
  )
  
  # PLOT
  ggplot(corr_df, aes(Var1, Var2, fill = Correlation)) +
    
    geom_tile(color = "white", size = 0.3) +
    
    # option valeurs (léger)
    geom_text(aes(label = ifelse(abs(Correlation) > 0.5,
                                 sprintf("%.2f", Correlation), "")),
              size = 3,
              color = "white",
              fontface = "bold") +
    
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    
    scale_x_discrete(labels = labels) +
    scale_y_discrete(labels = labels) +
    
    coord_fixed() +
    
    theme_minimal(base_size = 13) +
    
    labs(
      title = "Trait correlation matrix",
      subtitle = "Ordered by correlation with Model (GH marker)",
      fill = "r"
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, face = "bold"),
      axis.text.y.left = element_markdown(face = "bold"),
      
      axis.title = element_blank(),
      
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      
      legend.position = "right"
    )
}


# ==== 5. DIVERSITY PATTERNS ====

# - Comparison of local diversity patterns between above and belowground inventories
# - Correlations de diversité selon le degrée d’agregation

# ===== < DEPRECATED > =====



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



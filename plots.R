
# ===== plots.R =====
# Regroup plot functions by Subject

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr) # select, filter
library(ggplot2) # Plot
library(ggtext)
library(ggpubr) # afficher correlation sur graph
library(viridis) # turbo
library(patchwork)
library(akima)  # pour l'interpolation 2D
library(emmeans)
library(multcomp)
library(randomForest)


# ==== 1. TAXONOMY ====

# -- function : take two tables and make plot log-log --
plot_log_log <- function(Ab1, Ab2, marker){
  # JOIN
  Ab_join <- inner_join(Ab1, Ab2, by = "Family")

  # PLOT
  ggplot(Ab_join, aes(x = log(Abundance+1), y = log(Reads+1))) +
    theme_light() +
    geom_point(aes(color=tree)) +
    geom_smooth(method = "lm", color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_fixed() +
    xlim(0, 15) +
    ylim(0, 15) +
    stat_cor(method = "kendall",
             cor.coef.name = "tau",
             label.x = 1, label.y = 14) +
    annotate(
      "text",
      x = 1,
      y = 13,
      label = paste0("n = ", nrow(Ab_join)),
      hjust = 0,
      size = 4
    ) +
    labs(
      title = paste0(marker, " ~ Bota"),
      color = "trees (%)") +
    xlab("log(Basal Area + 1)")+
    ylab("log(Reads + 1)")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
    scale_color_gradient(
      low = "black",
      high = "blue",
      limits = c(0, 100)
    )
}

ab_rel_tab <- function(df_ab_rel){
  
  # Ordre des Sources
  df_ab_rel$Source <- factor(
    df_ab_rel$Source,
    levels = c("Ref","Bota","Sper01","ITS2","18S")
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
  
  # PLOT
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

fam_order_col <- function(df_list){
  
  ## Families share
  # All families
  all_fam <-  bind_rows(df_list) %>% 
    pull(Clade) %>%
    unique()
  
  # Ref families
  ref_fam <- df_list[5][[1]]$Clade
  
  # Bota families
  bota_fam <- df_list[4][[1]]$Clade
  
  # Genet families
  gen_fam <- bind_rows(df_list[c(1,2,3)]) %>% 
    pull(Clade) %>%
    unique()
  
  ## Growth form
  growth_form <- growth_form_family(all_fam)
  
  # Linéariser
  growth_form_long <- growth_form %>%
    pivot_longer(
      cols = -c(family),
      names_to = "growth_form",
      values_to = "Percent"
    )
  
  # Extraire l'ordre par tree percent
  tree_order <- growth_form_long %>%
    filter(growth_form == "tree") %>%
    arrange(desc(Percent)) %>%
    pull(family)
  
  # ajouter familles sans tree
  tree_order <- c(
    tree_order,
    setdiff(unique(growth_form_long$family), tree_order),
    "<2%","Indet."
  )
  
  ## FAMILIES REFERENCE
  all_fam_df <- data.frame(
    family = all_fam,
    Clade = all_fam,
    in_bota = all_fam %in% bota_fam & !all_fam %in% gen_fam,
    in_genet = all_fam %in% gen_fam & !all_fam %in% bota_fam,
    in_both = all_fam %in% bota_fam & all_fam %in% gen_fam
    ) %>% 
    mutate(
      in_ref = all_fam %in% ref_fam & !in_bota & !in_genet & !in_both
    )%>%
    mutate( # Ranks
      color_rank = case_when(
        Clade %in% c("<2%","Indet.") ~ 1,
        in_bota ~ 2,
        in_both ~ 3,
        in_genet ~ 4,
        in_ref  ~ 5
      ),
      tree_rank = match(family, tree_order),
      tree_rank = ifelse(is.na(tree_rank), max(match(family, tree_order), na.rm = TRUE) + 1, tree_rank),
      final_rank = color_rank*1000 + tree_rank # Final Rank
    ) %>%
    mutate( # Family color
      Label = case_when(
        Clade %in% c("<2%","Indet.") ~ paste0("<span style='color:black;'>", Clade, "</span>"),
        in_bota                      ~ paste0("<span style='color:forestgreen;'>", toupper(substr(Clade, 1, 3)), "</span>"),
        in_both                      ~ paste0("<span style='color:blue;'>", toupper(substr(Clade, 1, 3)), "</span>"),
        in_genet                     ~ paste0("<span style='color:red;'>", toupper(substr(Clade, 1, 3)), "</span>"),
        in_ref                       ~ paste0("<span style='color:brown;'>", toupper(substr(Clade, 1, 3)), "</span>"))
    ) %>%
    mutate( # fill color
      Color = case_when(
        Clade %in% c("<2%","Indet.") ~ "grey50",
        in_bota                      ~ "seagreen3",
        in_both                      ~ "steelblue2",
        in_genet                     ~ "indianred2",
        in_ref                       ~ "orange")
      )
  
  return(all_fam_df)
}

ab_rel_plot <- function(df_plot, treshold_label = 10){
  
  ### ORDER AND COLOR
  ## FAMILIES
  # Family labels order
  df_plot <- df_plot %>%
    mutate( 
      Clade = factor(
        Clade,
        levels = df_plot %>%
          distinct(Clade, final_rank) %>%
          arrange(final_rank) %>%
          pull(Clade)
      )
    )
  
  # Family labels color
  labels_vec <- setNames(df_plot$Label, df_plot$Clade)
  
  ## SOURCES
  # Source labels order
  df_plot <- df_plot %>%
    mutate(Source = factor(Source, levels = c("REF","Bota","Sper01","ITS2","18S")))
  
  # Source labels color
  labels_source <- c(
    "REF" = "<span style='color:orange;'> REF </span>",
    "Bota" = "<span style='color:forestgreen;'> Bota </span>",
    "Sper01" = "<span style='color:red;'> Sper01 </span>",
    "ITS2" = "<span style='color:red;'> ITS2 </span>",
    "18S" = "<span style='color:red;'> 18S </span>"
  )
  
  # Labels position
  df_plot <- df_plot %>%
    arrange(Source, desc(final_rank)) %>%
    group_by(Source) %>%
    mutate(
      ymax = cumsum(Metric),
      ymin = ymax - Metric,
      ymid = (ymin + ymax)/2
    ) %>%
    ungroup()
  
  # Fill color
  color_vec <- df_plot %>%
    distinct(Clade, Color) %>%
    deframe()
  
  # PLOT
  ggplot(df_plot, aes(x = Source, y = Metric, fill = Clade)) + 
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    
    ggtext::geom_richtext(
      data = df_plot %>% filter(Metric > treshold_label),
      aes(x = Source, y = ymid, label = Label),
      fill = "white",
      label.size = 0,
      fontface = "bold.italic",
      size = 3
    ) +
    
    scale_fill_manual(values = color_vec) +
    scale_x_discrete(labels = labels_source) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      x = "",
      y = "Proportion (%)",
      fill = "Clade"
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

growth_form_plot <- function(df_plot){
  
  ## SOURCES
  # Source labels order
  df_plot <- df_plot %>%
    mutate(Source = factor(Source, levels = rev(c("REF","Bota","Sper01","ITS2","18S"))))
  
  # Linéariser le dataframe
  growth_form_long <- df_plot %>%
    pivot_longer(
      cols = c(terrestrial,  tree, treelet, liana, shrub,  vine, epilithic, parasite, hemi.epiphyte),
      names_to = "growth_form",
      values_to = "Percent"
    )
  # Check 100%
  # growth_form_long %>% group_by(Source) %>% summarise(Percent = sum(Percent))
  
  # Source labels color
  labels_source <- c(
    "REF" = "<span style='color:orange;'> REF </span>",
    "Bota" = "<span style='color:forestgreen;'> Bota </span>",
    "Sper01" = "<span style='color:red;'> Sper01 </span>",
    "ITS2" = "<span style='color:red;'> ITS2 </span>",
    "18S" = "<span style='color:red;'> 18S </span>"
  )
  
  ## Growth Order
  # Import
  tab_ref <- tab_all_growth()
  
  # Apply
  growth_form_long$growth_form <- factor(
    growth_form_long$growth_form,
    levels = rev(tab_ref$type)
  )
  
  # Labels
  labels_growth <- paste0(
    tab_ref$type,
    " (",
    round(tab_ref$V1, 0),
    "%)"
  )
  names(labels_growth) <- tab_ref$type
  
  # PLOT
  ggplot(growth_form_long, aes(x = Percent, y = Source, fill = growth_form)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    scale_y_discrete(
      labels = labels_source,
      guide = guide_axis(n.dodge = 1)
    ) +
    scale_fill_viridis(
      discrete = TRUE,
      option = "C",
      labels = labels_growth,
      guide = guide_legend(reverse = TRUE),
      begin = 0.1,
      end = 0.9
    )+
    labs(
      x = "Proportion (%)",
      y = "",
      fill = "Growth Form"
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

plot_bias_growth <- function(df_bias, ref_tab){
  
  # Order factors
  ## inventory
  df_bias$inventory <- factor(
    df_bias$inventory,
    levels = c("18S","ITS2","Sper01","BOTA")
  )
  ## growth form
  df_bias$growth <- factor(
    df_bias$growth,
    levels = ref_tab$type
  )
  
  # Plots
  df_bias %>%
    ggplot(aes(x = growth, y = inventory, fill = bias)) +
    geom_tile() +
    geom_text(aes(label = round(bias, 2)), size = 4) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 1,
      limits = c(0, max(df_bias$bias, na.rm = TRUE))
    ) +
    theme_minimal() +
    theme(
      axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, face = "bold"),
      axis.text.y.left = element_markdown(face = "bold"),
    ) +
    labs(x = "Growth form", y = "Inventory", fill = "Ratio (inv/ref)")
}



# ===== 2. SPATIAL PATTERNS =====

plot_correlogramme <- function(tab_rad, Marker){
  
  # As Factors
  tab_rad$b.Radius <- as.factor(tab_rad$b.Radius)
  tab_rad$g.Marker <- as.factor(tab_rad$g.Marker)
  tab_rad$Level <- as.factor(tab_rad$Level)
  
  # Median and quartiles
  tab <- tab_rad %>%
    filter(g.Marker == Marker)
  
  # PLOT
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

corr_tab_plot <- function(config_tab){
  
  # SET ORDER
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("T+","T-","MCC","Count","EXP")
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
    levels = c("Sper01","ITS2","18S")
  )
  
  # clean data
  config_tab_clean <- droplevels(na.omit(config_tab))
  
  # PLOT
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

corr_violin_plot <- function(config_tab){
  
  config_tab <- config_tab %>%
    dplyr::select(g.Marker, Level, Metric_corr, Name, Correlation)
  
  # ordre des colonnes (Metric)
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("EXP","HYP","Count","T+","T-","MCC")
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
  
  # PLOT
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

morpho_corr_matrix <- function(level, marker, model){
  
  # Morpho
  morpho <- import_plants_morpho() %>%
    mutate(Species = gsub("_", " ", Species))
  
  # Correlations
  model_corr <- tab_rar %>%
    filter(Level == level,
           g.Marker == marker,
           Metric_corr == model) %>%
    dplyr::select(Name, Correlation) %>%
    rename(Species = Name, EXP = Correlation)
  
  # Join
  df <- inner_join(model_corr, morpho, by = level) %>%
    dplyr::select(-Species, -Genus, -Family, -order, -subclass2, -subclass1)
  
  # CORRELATION
  corr_mat <- cor(df, use = "pairwise.complete.obs")
  
  # Ordre basé sur Model
  order_vars <- names(sort(corr_mat[model, ], decreasing = TRUE))
  order_vars <- c(model, setdiff(order_vars, model))
  
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
    black = c(model)
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
      fill = expression("Pearson's " * italic(r))
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      #axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, face = "bold"),
      axis.text.y.left = element_markdown(face = "bold"),
      
      axis.title = element_blank(),
      
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      
      legend.position = "right",
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}


# ==== 3. DIVERSITY PATTERNS ====

# Random forest function
rand_for_IncMSE <- function(data_pixel){
  
  set.seed(123) # For reproductibility
  
  rf_model <- randomForest(
    Metric ~ .,              # toutes les variables explicatives
    data = data_pixel %>% dplyr::select(-x, -y, -sample_loc),  # enlever coordonnées si non pertinentes
    ntree = 500,
    importance = TRUE
  )
  
  # # importance des variables
  # importance(rf_model)
  # varImpPlot(rf_model)
  
  # Plot
  var_imp <- as.data.frame(importance(rf_model))
  var_imp$Variable <- rownames(var_imp)
  
  p <- ggplot(var_imp, aes(x = reorder(Variable, `%IncMSE`), y = `%IncMSE`)) +
    geom_col(fill = "grey") +
    coord_flip() +
    theme_light() +
    #scale_y_continuous(limits = c(0, 20)) + # Homogène scaling
    theme(
      plot.title = element_text(hjust = 0,size = 14, face="bold"),
      axis.text.y.left = element_text(face="bold")
    )
  
  return(list(p, rf_model, data_pixel))
}

# --- Accumulated Local Effects ---
plot_res_by_features <- function(res, features){
  
  library(iml)
  
  plot_list <- lapply(features, function(f){
    predictor <- Predictor$new(res[[2]], data = res[[3]], y = res[[3]]$Metric)
    ale <- FeatureEffect$new(predictor, feature = f, method = "ale")
    plot(ale)
  })
  
  return(plot_list)
}

# ==== < Basic > ====

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

plot_grid <- function(grid, stems = FALSE) {
  
  # Configuration
  df <- grid$df
  C <- grid$config
  
  # Création du ggplot de base
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_tile(aes(fill = Metric)) +
    coord_fixed() +
    theme_minimal() +
    labs(fill = NULL) +
    scale_fill_viridis_c(
      option = "turbo",
      breaks = c(min(df$Metric, na.rm = TRUE),
                 max(df$Metric, na.rm = TRUE)),
      labels = function(x) {
        toupper(format(x, scientific = TRUE, digits = 1))
      },
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,   # centre horizontalement
        face = "bold"  # texte en gras
      ),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",  # légende sous le plot
      legend.margin = ggplot2::margin(t=-10)
    )
  
  # Ajouter les stems si demandé
  
  if(stems) {
    bota_focus <- import_bota(list("Year" = "2012")) %>%
      dplyr::filter(.data[[C$Level]] == C$Name)
    p <- p +
      geom_point(
        data = bota_focus,
        aes(x = Xfield, y = Yfield, size = CircCorr/10),
        shape = 21,   # cercle avec contour
        fill = NA,    # pas de remplissage
        colour = "white",
        stroke = 0.5
      ) +
      guides(size = "none") +
      geom_point(
        data = bota_focus,
        aes(x = Xfield, y = Yfield),
        size = 4,
        shape = 3,   # cross
        colour = "white"
      )
  }
  
  return(p)
}

# Give a list of plot from a list of grid
grids_to_plots <- function(grid_list, stems = FALSE){
  
  plot_list <- lapply(seq_along(grid_list), function(i){
    plot_grid(grid_list[[i]], stems) +
      labs(title = names(grid_list)[i])
  })
  
  return(plot_list)
}

# ==== < Patch > ====

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
    Marker = "Sper01",
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

mono_turbo_species <- function(grid, corr, factor = 2) {
  
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
      plot.title = ggplot2::element_text(
        hjust = 0.5,   # centré
        size = 12,
        face = "bold.italic"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,   # centré
        size = 11,
        face = "bold"
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
    # ajouter nom centré en bas
    labs(
      title = grid$config$Name,
      caption = paste0("Pearson's r = ", corr)
    )
}


patch_turbo_species_ex <- function(N_patch){
  load("~/Bureau/ADNe_Nourague_2026/ADNe_Nourague_2026_files/.df_corr_tot.RData",
       envir = .GlobalEnv)
  
  common_sp <- df_corr_tot %>%
    dplyr::filter(
      Level %in% c("Species"),
      g.Marker %in% c("Sper01"),
      Metric_corr %in% c("EXP")
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
      "Marker" = "Sper01",
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
    p <- mono_turbo_species(genet_grid, corr)
    
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
    # p <- p +
    #   annotate(
    #     "text",
    #     x = (min(genet_grid$df$x)+max(genet_grid$df$x))/2,
    #     y = min(genet_grid$df$y - 60),
    #     label = paste0("Pearson's r = ", corr),
    #     fontface = "bold",
    #     hjust = 0.5,
    #     vjust = 1,
    #     size = 4,
    #     color = "black"
    #   )
  })
  
  # assembler les plots en patchwork
  wrap_plots(common_sp$plots, ncol = N_patch)
  # plot_annotation(
  #   title = Year,
  #   theme = theme(
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  #   )
  # )
}

patch_turbo_species_years <- function(Year1, Year2){
  
  # Import Genet species
  genet_sp <- give_genet_taxa(Rar = T, Level = "Species", Marker = "Sper01")
  
  # Import bota species
  bota1 <- import_bota(list("Year" = Year1))
  bota2 <- import_bota(list("Year" = Year2))
  
  # Keep unique species
  bota_sp1 <- unique(bota1$Species)
  bota_sp2 <- unique(bota2$Species)
  
  # Common species with genet
  common_sp <- Reduce(intersect, list(genet_sp, bota_sp1, bota_sp2)) %>%
    sort()
  
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
      "Marker" = "Sper01",
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
  #     subtitle = "Au niveau specifique pour le marqueur Sper01",
  #     theme = theme(
  #       plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #       plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12)
  #     )
  #   ) 
  
}

# ==== < Unused but useful > ====

# Whittaker plot
plot_Whittaker_bota <- function(df, Level){
  Ab_tib <- rank_ab_tib_bota(df, Level)
  Ab_vec <- rank_ab_vec_bota(Ab_tib, Level)
  
  plot(Ab_vec, Distribution = "lnorm", main = "Bota")
  param <- RAClnorm(Ab_vec)
  text(0.8*max(Ab_tib$Rank),0.5*max(Ab_tib$Abundance),paste0("mu = ", round(param$mu,digits = 3), "\n sigma = ", round(param$sigma,digits = 3)))
}
# Whittaker plot
plot_Whittaker_genet <- function(df, Level){
  Ab_tib <- rank_ab_tib_genet(df, Level)
  Ab_vec <- rank_ab_vec_genet(Ab_tib, Level)
  
  plot(Ab_vec, Distribution = "lnorm", ylab = "Reads")
  param <- RAClnorm(Ab_vec)
  text(0.8*max(Ab_tib$Rank),0.5*max(Ab_tib$Reads),paste0("mu = ", round(param$mu,digits = 3), "\n sigma = ", round(param$sigma,digits = 3)))
}

morpho_corr_bar <- function(marker = "Sper01", level_filter = "Species") {
  
  # Import Morpho
  morpho <- import_plants_morpho() %>%
    mutate(Species = gsub("_", " ", Species))
  
  # Read Model
  model <- tab_rar %>%
    dplyr::filter(
      Level == level_filter,
      g.Marker == marker,
      Metric_corr == "HYP"
    ) 
  %>%
    dplyr::select(Name, Correlation) %>%
    rename(Species = Name, Model = Correlation)
  
  # Inner Join
  model_morpho_df <- inner_join(model, morpho, by = "Species") %>%
    dplyr::select(
      -Species, -Genus, -Family, -order, -subclass2, -subclass1
    )
  
  # Correlation
  corr_vec <- cor(model_morpho_df, use = "pairwise.complete.obs")[, "Model"]
  
  corr_df <- data.frame(
    Variable = names(corr_vec),
    Correlation = corr_vec
  ) %>%
    dplyr::filter(Variable != "Model") %>%
    arrange(Correlation)
  
  # Colors
  colors <- list(
    forestgreen = c("Chlorophyll_content", "Thickness", "Toughness",
                    "Leaf_Area", "SLA", "C", "N", "13C", "Ca", "P", "K"),
    orange = c("Trunk_bark_thickness", "Sapwood_WSG"),
    brown = c("WSG", "Diameter", "SRL", "Tissue_density", "SRTA", "branchiness")
  )
  
  param_colors <- setNames(
    rep(names(colors), lengths(colors)),
    unlist(colors)
  )
  
  labels_source <- setNames(
    paste0("<span style='color:", param_colors, ";'> ", names(param_colors), " </span>"),
    names(param_colors)
  )
  
  # Plot
  ggplot(corr_df, aes(
    x = reorder(Variable, Correlation),
    y = Correlation,
    fill = Correlation
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Corrélation des traits morphologiques avec le modèle",
      subtitle = paste0("niveau spécifique pour le marker ", marker),
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

# ===== < DEPRECATED > =====



multi_filter_flow <- function(){
  
  # Init
  markers <- c("Sper01","18S","ITS2")
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
      mutate(Marker = factor(Marker, levels = c("Sper01", "18S", "ITS2")))
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



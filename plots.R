
# ===== plots.R =====
# Builds the figures of the report, grouped by subject.
# Each function returns a ggplot object (or a patchwork of them).

# ===== Libraries =====
library(tidyverse)     # Data manipulation and ggplot2
library(ggtext)        # Colored labels (element_markdown)
library(ggpubr)        # Correlation displayed on a plot (stat_cor)
library(viridis)       # Color palettes
library(patchwork)     # Arrange plots
library(akima)         # 2D interpolation
library(emmeans)       # Adjusted means after the ANOVA
library(multcomp)      # Letters of significance (cld)
library(randomForest)  # Random forests
library(ggVennDiagram) # Venn diagrams
library(UpSetR)        # UpSet diagrams

# Morphological traits by organ, the name of each group is its color in plots
TRAIT_COLORS <- list(
  forestgreen = c("Chlorophyll_content", "Thickness", "Toughness",
                  "Leaf_Area", "SLA", "C", "N", "13C", "Ca", "P", "K"), # Leaves
  orange = c("Trunk_bark_thickness", "Sapwood_WSG"),                    # Stems
  brown = c("WSG", "Diameter", "SRL", "Tissue_density", "SRTA", "branchiness") # Roots
)

# ==== 1. TAXONOMY ====

#' Plot the reads against the basal area of the families (Figure 4)
#'
#' Each panel gives the Kendall tau, the p-value of the interaction with the
#' percentage of trees in a linear model, and the number of families.
#'
#' @param df_H2 Data frame with Family, Basal_Area, tree (percentage of tree
#'   species in the family), Reads and Marker.
#' @return A ggplot with one panel per marker.
plot_log_log <- function(df_H2){

  # Title labels
  lab_title <- list(
    "Sper01" = "Sper01 ~ Bota",
    "ITS2" = "ITS2 ~ Bota"
  )
  
  # Factorize
  df_H2$Marker <- factor(
    df_H2$Marker,
    levels = c("Sper01","ITS2")
  )
  
  # compute n per marker
  n_df <- df_H2 %>%
    dplyr::group_by(Marker) %>%
    dplyr::summarise(n = n(), .groups = "drop")
  
  # Linear model with the tree interaction, for each marker
  model_df <- df_H2 %>%
    group_by(Marker) %>%
    tidyr::nest() %>%
    mutate(
      model = purrr::map(data, ~ lm(
        log(Reads + 1) ~ log(Basal_Area + 1) * tree,
        data = .x
      ))
    )
  # p-value of the interaction
  p_table <- model_df %>%
    mutate(
      anova = purrr::map(model, ~ car::Anova(.x, type = 2)),
      p_interaction = purrr::map_dbl(anova, ~ .x["log(Basal_Area + 1):tree", "Pr(>F)"])
    ) %>%
    dplyr::select(Marker, p_interaction)
  
  # PLOT
  ggplot(df_H2, aes(x = log(Basal_Area+1), y = log(Reads+1))) +
    # facet_wrap
    facet_wrap(~Marker, scales = "fixed", strip.position = "top", labeller = labeller(Marker = lab_title))+
    geom_point(aes(color=tree)) +
    geom_smooth(method = "lm", color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_fixed() +
    xlim(0, 15) +
    ylim(0, 15) +
    # Stat
    stat_cor(
      method = "kendall",
      cor.coef.name = "tau",
      label.x = 1,
      label.y = 14
    ) +
    geom_text(
      data = p_table,
      aes(x = 1, y = 12.5, label = paste0("p(tree) = ", signif(p_interaction, 2))),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4
    )+
    # n per facet
    geom_text(
      data = n_df,
      aes(x = 1, y = 11, label = paste0("n = ", n)),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4
    ) +
    labs(color = "trees (%)") +
    xlab("log(Basal Area + 1)")+
    ylab("log(Reads + 1)")+
    scale_color_gradient(
      low = "black",
      high = "blue",
      limits = c(0, 100)
    ) +
    theme_bw() +
    theme(strip.text.x = element_text(face = "bold",size=14, colour = "black"))
}

#' Plot the table of the relative richness by family (Figure S2)
#'
#' @param df_ab_rel Data frame with Clade, Metric (percentage) and Source.
#' @return A ggplot table, percentages above 2% are in bold.
ab_rel_tab <- function(df_ab_rel){
  
  # Order of the sources
  df_ab_rel$Source <- factor(
    df_ab_rel$Source,
    levels = c("Ref","Bota","Sper01","ITS2","18S")
  )
  
  # Order of the clades: by decreasing share in Bota
  order_bota <- df_ab_rel %>%
    filter(Source == "Bota") %>%
    arrange(desc(Metric)) %>%
    distinct(Clade) %>%
    pull(Clade)
  
  # Other clades (not in Bota)
  other_clade <- setdiff(unique(df_ab_rel$Clade), order_bota)
  
  # Combine both
  order_clade <- c(order_bota, other_clade)
  
  # Apply the order
  df_ab_rel$Clade <- factor(df_ab_rel$Clade, levels = rev(order_clade))
  
  # Clean data
  tab_clean <- df_ab_rel %>%
    dplyr::filter(!is.na(Metric)) %>%
    droplevels()
  
  # Y labels, bold above 2%
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

#' Compute the display order and the colors of the families (Figure 3)
#'
#' Families are grouped by the surveys that detect them, then sorted by their
#' percentage of tree species.
#'
#' @param df_list List of 5 data frames in the order of family_richness_list().
#' @return A data frame with one row per family: Clade, the presence flags
#'   (in_bota, in_genet, in_both, in_ref), final_rank, Label and Color.
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
  
  # Long format
  growth_form_long <- growth_form %>%
    pivot_longer(
      cols = -c(family),
      names_to = "growth_form",
      values_to = "Percent"
    )
  
  # Order by percentage of trees
  tree_order <- growth_form_long %>%
    filter(growth_form == "tree") %>%
    arrange(desc(Percent)) %>%
    pull(family)
  
  # Add the families without trees
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

#' Plot the relative richness of the families in each survey (Figure 3a)
#'
#' @param df_plot Data frame with Clade, Metric, Source, and the columns of
#'   fam_order_col().
#' @param treshold_label Percentage above which a family is labelled.
#' @return A ggplot of stacked bars.
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

#' Plot the share of the growth forms in each survey (Figure 3b)
#'
#' @param df_plot Result of growth_form_source().
#' @return A ggplot of stacked bars.
growth_form_plot <- function(df_plot){
  
  ## SOURCES
  # Source labels order
  df_plot <- df_plot %>%
    mutate(Source = factor(Source, levels = rev(c("REF","Bota","Sper01","ITS2","18S"))))
  
  # Long format
  growth_form_long <- df_plot %>%
    pivot_longer(
      cols = c(terrestrial,  tree, treelet, liana, shrub,  vine, epilithic, parasite, hemi.epiphyte),
      names_to = "growth_form",
      values_to = "Percent"
    )
  
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

# ===== 2. SPATIAL PATTERNS AND TRAITS =====

#' Plot the correlation of the diversity patterns by radius (Figure 9)
#'
#' @param tab_rad Result of make_cor_table() with a column b.Radius.
#' @return A ggplot with one panel per marker.
plot_correlogramme <- function(tab_rad){
  
  # As Factors
  tab_rad$b.Radius <- as.factor(tab_rad$b.Radius)
  tab_rad$g.Marker <- factor(
    tab_rad$g.Marker, 
    levels = c("Sper01","ITS2","18S"))
  tab_rad$Level <- as.factor(tab_rad$Level)
  
  # Significance of the correlations
  tab <- tab_rad %>%
    filter(g.Marker %in% c("Sper01","ITS2")) %>%
    mutate(
      p.value = ifelse(p.value <= 0.05, "p<=0.05", "p>0.05")
    )
  
  # PLOT
  ggplot(tab, aes(x = b.Radius, y = Correlation, color = Level, group=Level)) +
    facet_wrap(~g.Marker, scales = "fixed", strip.position = "top")+
    geom_line() +
    geom_point(aes(shape = p.value), size = 2.5) +
    scale_shape_manual(
      values = c(
        `p<=0.05` = 16, # Filled circle
        `p>0.05` = 1    # Empty circle
      )
    )+
    labs(x = "Radius (m)", y = "Spearmans' r") +
    theme_bw() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    theme(strip.text.x = element_text(face = "bold",size=14))
}

#' Plot the correlations between botanical and eDNA patterns (Figure 7)
#'
#' Letters come from pairwise comparisons of the levels after a three-way ANOVA
#' (level, method, marker), with a threshold of 0.1.
#'
#' @param config_tab Result of make_cor_table().
#' @return A ggplot of violins with one panel per marker and level.
corr_violin_plot <- function(config_tab){
  
  config_tab <- config_tab %>%
    dplyr::select(g.Marker, Level, Metric_corr, Name, Correlation)
  
  # Metric
  config_tab$Metric_corr <- factor(
    config_tab$Metric_corr,
    levels = c("EXP","HYP","Count","I.O")
  )
  
  # Level
  config_tab$Level <- factor(
    config_tab$Level,
    levels = c("Family","Genus","Species")
  )
  
  # Marker
  config_tab$g.Marker <- factor(
    config_tab$g.Marker,
    levels = c("Sper01","ITS2")
  )
  
  # clean data
  config_tab_clean <- config_tab %>%
    filter(!is.na(Correlation)) %>%
    droplevels()
  
  # Three-way ANOVA
  res_aov3 <- aov(Correlation ~ Level * Metric_corr * g.Marker, 
                  data = config_tab_clean)
  
  # Adjusted means by level
  emm3 <- emmeans(res_aov3, ~ Level | Metric_corr * g.Marker)
  
  # Letters of significance
  cld3 <- cld(emm3, Letters = letters, alpha = 0.1) # p-value 0.1
  cld_df3 <- as.data.frame(cld3)
  
  # Position of the letters
  cld_df3 <- cld_df3 %>%
    left_join(
      config_tab %>%
        group_by(Metric_corr, g.Marker, Level) %>%
        summarise(
          pos = min(Correlation, na.rm = TRUE)
        ),
      by = c("Metric_corr", "g.Marker", "Level")
    )
  
  # Number of points
  n_df <- config_tab_clean %>%
    group_by(Metric_corr, g.Marker, Level) %>%
    summarise(n = n(), .groups = "drop")
  
  n_df <- n_df %>%
    mutate(pos = 0.85)
  
  # PLOT
  ggplot(config_tab_clean,
         aes(x = Correlation,
             y = Metric_corr)) +
    
    # violins
    geom_violin(
      width = 0.8,
      fill = "grey80",
      color = "black"
    ) +
    
    # points
    geom_jitter(
      width = 0.15,
      size = 1,
      alpha = 0.5
    ) +
    
    # facets
    facet_grid(g.Marker ~ Level,
               scales = "free_x",
               space = "free_x",
               switch = "y",
               drop = TRUE) +
    
    geom_text(data = cld_df3,
              aes(x = pos - 0.05,
                  y = Metric_corr,
                  label = .group),
              inherit.aes = FALSE,
              size = 5,
              col = "red") +
    geom_text(
      data = n_df,
      aes(x = pos,
          y = Metric_corr,
          label = paste0("n=", n)),
      inherit.aes = FALSE,
      size = 4
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 1, linetype = "dashed", color = "black")+
    
    theme_minimal(base_size = 11) +
    
    labs(
      x = "Correlation coefficient"
    ) + 
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(-0.15, 1.1)
    )+
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x.bottom = element_text(size=13),
      axis.line = element_line(color = "gray"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),

      # strips
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(face = "bold",size=14),
      strip.text.y.left = element_text(angle = 0, face = "bold",size=14),

      # axes
      axis.text.x.top = element_text(face = "bold"),
      axis.text.x.bottom = element_text(size = 13),
      axis.text.y = element_text(size = 14)
    )
}

#' Correlate the eDNA-botanical correlations with the morphological traits
#'
#' @param df_corr Result of make_cor_table().
#' @param marker,level,metric Selection of the rows of df_corr.
#' @param corr Method of cor.test(), for example "pearson".
#' @return A data frame with Variable (trait), Correlation, p.value, signif
#'   (stars) and Level.
morpho_corr_df <- function(df_corr, marker = "Sper01", level = "Species", metric = "HYP", corr) {
  
  morpho <- import_plants_morpho() %>%
    mutate(Species = gsub("_", " ", Species))
  
  morpho_agg <- aggregate_morpho(morpho, level, median)
  
  model <- df_corr %>%
    dplyr::filter(
      Level == level,
      g.Marker == marker,
      Metric_corr == metric
    ) %>%
    dplyr::select(Name, Correlation) %>%
    rename(!!level := Name, Model = Correlation)
  
  model_morpho_df <- inner_join(model, morpho_agg, by = level) %>%
    dplyr::select(-!!level, -n)
  
  res_list <- lapply(names(model_morpho_df), function(var) {
    
    if (var == "Model") return(NULL)
    
    test <- cor.test(
      model_morpho_df[[var]],
      model_morpho_df$Model,
      method = corr,
      use = "pairwise.complete.obs"
    )
    
    data.frame(
      Variable = var,
      Correlation = test$estimate,
      p.value = test$p.value
    )
  })
  
  res_df <- bind_rows(res_list)
  
  res_df <- res_df %>%
    mutate(
      signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Level = level
    )
  
  return(res_df)
}

#' Plot the correlations with the morphological traits by level (Figure 8)
#'
#' @param df_corr Result of make_cor_table().
#' @param Marker,Metric Selection of the rows of df_corr.
#' @param Corr Method of cor.test(), for example "pearson".
#' @return A ggplot of bars with one panel per taxonomic level.
morpho_corr_bar_strip <- function(df_corr, Marker = "Sper01", Metric = "HYP", Corr){
  
  levels <- c("Family", "Genus", "Species")
  
  corr_all <- map_dfr(levels, ~ morpho_corr_df(df_corr, level = .x, marker = Marker, metric = Metric, corr = Corr))
  
  # Colors
  labels_source <- trait_labels()
  
  # Reorder
  family_order <- corr_all %>%
    filter(Level == "Family") %>%
    arrange(Correlation) %>%
    pull(Variable)
  
  corr_all <- corr_all %>%
    mutate(
      Variable = factor(Variable, levels = family_order)
    )
  
  # Labels
  n_df <- df_corr %>%
    group_by(Metric_corr, g.Marker, Level) %>%
    summarise(n = n(), .groups = "drop") %>%
    dplyr::filter(g.Marker == Marker, Metric_corr == Metric) %>%
    mutate(label = paste0(Level, " (n=", n, ")")) %>%
    dplyr::select(Level, label) %>%
    deframe()
  
  # Plot
  ggplot(corr_all, aes(
    x = Correlation,
    y = Variable,
    fill = Correlation
  )) +
    coord_cartesian(xlim = c(-1, 1))+
    geom_col() +
    scale_y_discrete(labels = labels_source) +
    # facet_wrap
    facet_wrap(~Level, scales = "fixed", strip.position = "top", labeller = labeller(Level = n_df))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.3),

      axis.title.y = element_blank(),
      axis.title.x.bottom = element_text(size=13),
      axis.line = element_line(color = "gray"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      
      # strips
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 14),
      
      # axes
      axis.text.x.bottom = element_text(size = 11),
      axis.text.y.left = ggtext::element_markdown(hjust = 1, face = "bold",size = 11),
      
      title = element_text(face = "bold", size = 16),
    ) +
    labs(
      x = Corr,
      fill = "r"
    ) +
    geom_text(
      aes(label = signif),
      size = 4,
      fontface = "bold",
      color = "black",
      nudge_x = 0.05
    )
}

#' Build the colored labels of the morphological traits
#'
#' @return A named vector of HTML labels, for ggtext::element_markdown().
trait_labels <- function(){
  
  param_colors <- setNames(
    rep(names(TRAIT_COLORS), lengths(TRAIT_COLORS)),
    unlist(TRAIT_COLORS)
  )
  
  setNames(
    paste0("<span style='color:", param_colors, ";'> ", names(param_colors), " </span>"),
    names(param_colors)
  )
}

#' Plot the Venn diagrams of the taxa shared between surveys (Figure 2)
#'
#' @param list_by_level List of taxa: list_by_level[[Level]][[Survey]].
#' @return A patchwork with one Venn diagram per taxonomic level.
plot_venn_levels <- function(list_by_level){
  
  plot_list <- lapply(names(list_by_level), function(lvl){
    
    col <- switch(lvl,
                  "Family"  = "darkorange",
                  "Genus"   = "steelblue",
                  "Species" = "forestgreen")
    
    ggVennDiagram(list_by_level[[lvl]]) +
      scale_fill_gradient(low = "white", high = col) +
      ggtitle(lvl) +
      guides(fill = guide_legend(
        override.aes = list(size = 10)
      )) +
      theme(
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "none"
      )
  })
  
  wrap_plots(plot_list, ncol = 3)
}

#' Plot the UpSet diagram of the taxa shared between surveys (Figure S1)
#'
#' @param tax_list List of taxa of each survey, at one taxonomic level.
#' @param level Name of the level, used in the axis title.
#' @param size_max Upper limit of the axis of the set sizes.
#' @return An UpSetR plot.
plot_upset_level <- function(tax_list, level, size_max){
  
  data_upset <- fromList(tax_list)
  
  upset(data_upset,
        nsets = ncol(data_upset),
        order.by = "freq",
        show.numbers = "yes",
        set_size.show = TRUE,
        main.bar.color = "steelblue",
        sets.bar.color = "gray",
        mainbar.y.label = paste("Intersections -", level),
        sets.x.label = "Set Size",
        set_size.scale_max = size_max)
}

#' Plot the distribution of the morphological traits (Figure S3)
#'
#' The stars after the name of a trait give the significance of the Shapiro
#' normality test.
#'
#' @param morpho Traits of the species, see import_plants_morpho().
#' @return A ggplot with one histogram per trait.
plot_trait_distribution <- function(morpho){
  
  morpho_long <- morpho %>%
    dplyr::select(where(is.numeric)) %>%
    pivot_longer(where(is.numeric), values_to = "value", names_to = "var")
  
  # Shapiro test on each trait
  shapiro_results <- morpho_long %>%
    group_by(var) %>%
    summarise(
      p_value = shapiro.test(value)$p.value
    ) %>%
    mutate(
      signif = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ "")
    ) %>%
    mutate(label = paste0(var, signif)) %>%
    dplyr::select(var, label) %>%
    deframe()
  
  ggplot(morpho_long, aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ var,
               scales = "free",
               ncol = 4,
               labeller = labeller(var = shapiro_results)) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"))
}

#' Plot the upper triangle of a correlation matrix (Figures S4 and S5)
#'
#' @param corr_mat Square correlation matrix.
#' @param threshold The values are displayed when their absolute value is above.
#' @param y_labels Labels of the variables on the y axis, or waiver() for names.
#' @param fill_label Title of the color legend.
#' @return A ggplot.
plot_corr_triangle <- function(corr_mat, threshold, y_labels = waiver(),
                               fill_label = expression("Pearson's " * italic(r))){
  
  # Long format
  corr_df <- as.data.frame(as.table(corr_mat))
  colnames(corr_df) <- c("Var1", "Var2", "Correlation")
  
  # Triangle (upper only)
  corr_df <- corr_df %>%
    dplyr::filter(as.numeric(Var1) <= as.numeric(Var2))
  
  ggplot(corr_df, aes(Var1, Var2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(abs(Correlation) > threshold,
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
    scale_y_discrete(labels = y_labels) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    labs(
      fill = fill_label
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.y.left = element_markdown(face = "bold"),
      axis.title = element_blank(),
      legend.position = "right",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# ==== 3. DIVERSITY PATTERNS ====

#' Estimate the importance of the variables with repeated random forests
#'
#' @param data_pixel Data frame with Metric (response variable), x, y,
#'   sample_loc (not used by the model) and the explanatory variables.
#' @param nrep Number of random forests, each one has 500 trees.
#' @return A long data frame with Variable, replicate and importance (%IncMSE).
rand_for_IncMSE_rep <- function(data_pixel, nrep = 100){
  
  results <- replicate(nrep, {
    
    rf_model <- randomForest(
      Metric ~ .,
      data = data_pixel %>% dplyr::select(-x, -y, -sample_loc),
      ntree = 500,
      importance = TRUE
    )
    
    importance(rf_model)[, "%IncMSE"]
    
  }, simplify = "matrix")
  
  # Long format
  df_imp <- as.data.frame(results)
  df_imp$Variable <- rownames(df_imp)
  
  df_long <- tidyr::pivot_longer(df_imp,
                                 cols = -Variable,
                                 names_to = "replicate",
                                 values_to = "importance")
  
  return(df_long)
}

#' Plot the importance of the variables in the diversity patterns (Figure 11)
#'
#' Importances are normalized within each random forest. Letters come from
#' Dunn tests between patterns for each variable (BH correction, p < 0.05).
#'
#' @param df_plot Results of rand_for_IncMSE_rep(), with a column pattern.
#' @return A ggplot of boxplots with one panel per variable.
plot_RF_stat <- function(df_plot){
  
  # Ordering
  df_plot$Variable <- recode(df_plot$Variable, "CNratio" = "C:N")
  
  df_plot$pattern <- factor(
    df_plot$pattern,
    levels = c("Bota","ITS2","Sper01")
  )
  
  df_plot$Variable <- factor(
    df_plot$Variable,
    levels = c("Al","Ca","Fe","C:N","K","Na","light","wetness")
  )
  
  # Normalize the importance within each random forest
  df_plot <- df_plot %>%
    group_by(pattern, replicate) %>%
    mutate(
      importance = importance / sum(importance)
    ) %>%
    ungroup()
  
  # Tests and letters
  letters_df <- df_plot %>%
    group_by(Variable) %>%
    do({
      
      # Kruskal
      kw <- kruskal.test(importance ~ pattern, data = .)
      
      # Dunn post-hoc
      dunn <- FSA::dunnTest(importance ~ pattern, data = ., method = "bh")
      
      # Lettres
      cld <- rcompanion::cldList(
        P.adj ~ Comparison,
        data = dunn$res,
        threshold = 0.05
      )
      
      cld$Variable <- unique(.$Variable)
      cld
    }) %>%
    ungroup() %>%
    dplyr::rename(pattern = Group) %>%
    # cldList() removes the zeros from the names of the groups
    mutate(pattern = recode(pattern,
                            "Sper1" = "Sper01"))
  
  # Plot
  ggplot(df_plot, aes(x = pattern, y = importance)) +
    geom_boxplot() +
    facet_wrap(~Variable, scales = "fixed") +
    geom_text(
      data = letters_df,
      aes(x = pattern, y = Inf, label = MonoLetter),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 5,
      color = "red"
    ) +
    theme_bw() +
    theme(
      strip.text.x = element_text(face = "bold",size=12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=13),
      axis.text.x = element_text(size=11,face="bold")
      ) +
    labs(y="Relative importance (%)")
}

# ==== Maps ====

#' Plot a grid as a heatmap (Figure 6)
#'
#' @param grid Grid, see grid_read().
#' @param stems If TRUE, add the stems of the taxon of the grid (2012 inventory)
#'   as white circles proportional to their circumference.
#' @return A ggplot.
plot_grid <- function(grid, stems = FALSE) {
  
  # Configuration
  df <- grid$df
  C <- grid$config
  
  # Base plot
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
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t=-10),
      legend.text = element_text(size = 11)
    )
  
  # Add the stems if requested
  if(stems) {
    bota_focus <- import_bota(list("Year" = "2012")) %>%
      dplyr::filter(.data[[C$Level]] == C$Name)
    p <- p +
      geom_point(
        data = bota_focus,
        aes(x = Xfield, y = Yfield, size = CircCorr/10),
        shape = 21,   # Empty circle
        fill = NA,
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

#' Interpolate a grid on a finer grid with bicubic splines
#'
#' @param grid Grid, see grid_read().
#' @param factor Number of new cells per cell of the grid, in each direction.
#' @return A data frame with x, y and krig (interpolated value, at least 0).
interpolate_grid_2D <- function(grid, factor = 3) {
  
  # Configuration
  G <- grid
  df <- G$df
  
  # 1 Fine grid
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
    linear = FALSE, # Bicubic spline
    extrap = TRUE
  )
  
  # 3 Data frame for ggplot
  interp_df <- expand.grid(x = interp_res$x, y = interp_res$y)
  interp_df$krig <- as.vector(interp_res$z)
  
  # 4 Truncate to 0
  interp_df$krig <- pmax(as.vector(interp_res$z), 0)
  
  return(interp_df)
}

#' Plot an interpolated grid with axes (Figure 10)
#'
#' @param grid Grid, see grid_read().
#' @param factor See interpolate_grid_2D().
#' @return A ggplot.
plot_grid_interp_soft <- function(grid, factor = 1) {
  
  fine <- interpolate_grid_2D(grid, factor)
  
  ggplot(fine, aes(x = x, y = y, fill = krig)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = turbo(100)
    ) +
    coord_fixed() +
    theme_light()+
    theme(
      legend.margin = ggplot2::margin(t = -10),
      plot.title = element_text(hjust = 0,size = 14, face="bold")
    )
}

#' Plot an interpolated grid without axes, with a scale bar (Figure 5)
#'
#' @param grid Grid, see grid_read().
#' @param factor See interpolate_grid_2D().
#' @return A ggplot titled with the name of the taxon.
mono_turbo_species <- function(grid, factor = 1) {
  
  fine <- interpolate_grid_2D(grid, factor)
  
  # Limits to place the scale bar
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
    
    # Remove axes and ticks
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t=-10),
      legend.text = element_text(size = 11),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold.italic"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        size = 11,
        face = "bold"
      ),
      plot.margin = ggplot2::margin(0, 20, 0, 20)
    ) +
    
    # Horizontal legend without title
    guides(
      fill = guide_colorbar(
        title = NULL, 
        direction = "horizontal",
        barwidth = unit(4, "lines"),
        barheight = unit(0.5, "lines")
      )
    ) +
    
    # Scale bar
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
      x = x_min + 20 + 25,
      y = y_max - 0.08 * (y_max - y_min),
      label = "50 m",
      hjust = 0.5,
      vjust = 1,
      size = 3,
      color = "white"
    ) +
    # Name of the taxon
    labs(
      title = grid$config$Name
    )
}

#' Plot the eDNA of some taxa with their stems in 2012 and 2022 (Figure 5)
#'
#' @param df_corr Result of make_cor_table(), it gives the taxa sorted by
#'   increasing correlation.
#' @param marker,level,metric Selection of the rows of df_corr.
#' @param selection Ranks of the taxa to plot, in this sorted list.
#' @param n Number of columns.
#' @return A patchwork with one map per taxon.
patch_turbo_ex <- function(df_corr, marker, level, metric, selection, n){
  
  # Level and Marker selection
  common_tax <- df_corr %>%
    dplyr::filter(Level == level, g.Marker == marker, Metric_corr == metric) %>%
    arrange(Correlation) %>%
    pull(Name)
  
  # Apply indices selection
  common_tax <- common_tax[selection]
  
  plots <- lapply(common_tax,function(sp){
    
    # Import genet grid
    genet_grid <- grid_read(list(
      "Type" = "genet",
      "Rar" = TRUE,
      "Marker" = marker,
      "Level" = level,
      "Name" = sp,
      "Metric" = "Reads"
    ))
    
    # Bota year 1 grid
    bota_focus1 <- import_bota(list("Year" = 2012)) %>%
      filter(.data[[level]] == sp) %>%
      mutate(Year = factor(2012, levels = c(2012, 2022)))
    
    # Bota year 2 grid
    bota_focus2 <- import_bota(list("Year" = 2022)) %>%
      filter(.data[[level]] == sp) %>%
      mutate(Year = factor(2022, levels = c(2012, 2022)))
    
    # Heatmap
    p <- mono_turbo_species(genet_grid)
    
    # Year 1
    p <- p +
      geom_point(
        data = bota_focus1,
        aes(x = Xfield, y = Yfield, size = CircCorr/4),
        color = "red",
        shape = 21,
        stroke = 0.5,
        fill = NA
      )
    
    # Year 2
    p <- p +
      geom_point(
        data = bota_focus2,
        aes(x = Xfield, y = Yfield, size = CircCorr),
        color = "orange",
        shape = 21,
        stroke = 0.5,
        fill = NA
      )+
      guides(size = "none")
    
    p
  })
  
  # Patchwork
  wrap_plots(plots, ncol = n)
}


# ===== table.R =====
# Builds the tables behind the figures of the report: taxa lists, relative
# richness by family, growth forms and correlations between surveys.

# ===== Libraries =====
library(tidyverse) # Data manipulation

# ===== 1. TAXONOMY =====

# ---- Taxa lists ----

#' List the taxa of a survey at a taxonomic level
#'
#' @param Level "Family", "Genus" or "Species".
#' @param Marker "Bota" (2012 inventory), "REF" (reference checklist) or an
#'   eDNA marker: "Sper01", "ITS2", "18S".
#' @return A vector of taxa names.
get_tax_list <- function(Level, Marker){
  # Bota
  if(Marker == "Bota"){
    df <- import_bota(list(Year = 2012))
    return(
      df %>%
        pull(.data[[Level]]) %>%
        unique()
    )
  }
  # Genet
  if(Marker %in% c("Sper01", "ITS2", "18S")){
    df <- import_genet(
      list(
        Marker = Marker,
        Rar = TRUE
      )
    )
    return(
      rank_ab_tib_genet(df, tolower(Level)) %>%
        pull(.data[[tolower(Level)]])
    )
  }
  # Ref
  if(Marker == "REF"){
    return(
      import_plants_nouragues() %>%
        pull(tolower(Level)) %>%
        unique()
    )
  }

  # else
  NULL

}

#' Keep the taxa shared by the botanical survey and each other survey
#'
#' @param list_by_level List of taxa: list_by_level[[Level]][[Survey]], with
#'   a survey named "Bota".
#' @return The same structure without "Bota": for each level and survey, the
#'   taxa also found in the botanical survey.
filter_common_taxa <- function(list_by_level){
  list_filtered <- lapply(names(list_by_level), function(level) {

    lvl_list <- list_by_level[[level]]

    bota <- lvl_list$Bota

    markers <- setdiff(names(lvl_list), "Bota")

    res <- lapply(markers, function(m) {
      intersect(bota, lvl_list[[m]])
    })

    names(res) <- markers
    res
  })

  names(list_filtered) <- names(list_by_level)

  return(list_filtered)
}

#' List the taxa of a survey from the parameters of a grid
#'
#' @param config List with Type ("bota" or "genet"), Level, and the parameters
#'   to import the survey (Year or Marker and Rar).
#' @return A vector of taxa names, without NA.
give_level_names <- function(config){

  if(config$Type == "genet"){
    # Import genet
    genet <- import_genet(config)

    # names in the level
    l <- tolower(config$Level) # tolower
    names <- unique(genet$motus[[l]]) # unique
    names <- names[!is.na(names)] # Remove NA
  }
  else if(config$Type == "bota"){
    # Import Bota
    bota <- import_bota(config)

    # names in the level
    names <- unique(bota[[config$Level]])
  }
  else{
    stop(paste0("Unknown config$Type : ", config$Type))
  }

  return(names)
}

# ---- Abundance by taxon ----

#' Compute the basal area of each taxon of a botanical inventory
#'
#' @param df Botanical inventory, see import_bota().
#' @param Level Name of the column with the taxonomic level.
#' @return A tibble with the taxon and Basal_Area, sorted by decreasing area.
rank_ab_tib_bota_BA <- function(df, Level){
  df %>%
    group_by(.data[[Level]]) %>%
    summarize(Basal_Area = sum(pi*(CircCorr/(2*pi))^2)) %>%
    arrange(desc(Basal_Area))
}

#' Compute the number of reads of each taxon of a marker
#'
#' @param df eDNA data, see import_genet().
#' @param Level Name of the column of df$motus with the taxonomic level
#'   ("family", "genus" or "species"). MOTUs without assignment are removed.
#' @return A tibble with the taxon, Reads, Dominance (proportion of the reads)
#'   and Rank, sorted by decreasing number of reads.
rank_ab_tib_genet <- function(df, Level){

  Reads <- colSums(df$reads)
  as.data.frame(Reads) -> Reads
  Reads$id <- row.names(Reads)

  left_join(df$motus, Reads, by = "id") %>%
    dplyr::filter(!is.na(.data[[Level]])) %>%
    group_by(.data[[Level]]) %>%
    summarise(Reads = sum(Reads))%>%
    ungroup() %>%
    mutate(Dominance = Reads / sum(Reads)) %>%
    arrange(desc(Reads)) %>%
    mutate(Rank = row_number())
}

# ---- Relative richness ----

#' Compute the share of each taxon in the total of a survey
#'
#' The Metric of the config is summed over the plot for all the taxa, then for
#' each taxon. Taxa without assignment are grouped as "Indet.".
#'
#' @param config List of the parameters of a grid, without Name.
#' @return A data frame with Clade (taxon) and Metric (percentage).
summary_ab_rel <- function(config){

  # Config All
  config_all <- config
  config_all$Name <- "All"

  # Read grid All
  grid_all <- grid_read(config_all)
  total <- sum(grid_all$df$Metric, na.rm = TRUE)

  # All the names of the clade
  all_names <- c(give_level_names(config),"Indet.") %>% unique()

  results <- lapply(all_names, function(n) {

    # Config Name
    config_name <- config
    config_name$Name <- n

    # Read data
    grid <- grid_read(config_name)

    # Total
    total_name <- sum(grid$df$Metric, na.rm = TRUE)
    ab_rel <- (total_name*100) / total

    data.frame(
      Clade = n,
      Metric = ab_rel
    )
  })

  final_df <- do.call(rbind, results)
  return(final_df)
}

#' Compute the share of each taxon for several surveys
#'
#' @param config_table Data frame with the parameters of a survey in each row.
#' @return config_table with two more columns: Result (list of the tables of
#'   summary_ab_rel()) and Check_100 (sum of the percentages, to check).
test_ab_rel_sum <- function(config_table){
  # Calcul tabs
  config_table$Result <- lapply(seq_len(nrow(config_table)), function(i) {
    message("Iteration ", i)
    summary_ab_rel(as.list(config_table[i, ]))
  })

  # Check 100%
  config_table$Check_100 <- lapply(seq_len(nrow(config_table)), function(i) {
    sum(config_table$Result[[i]]$Metric, na.rm = TRUE)
  })

  return(config_table)
}

#' Gather the richness by family of each survey and of the reference
#'
#' @param df_H1 Result of test_ab_rel_sum() for Sper01, ITS2, 18S and Bota, in
#'   this order.
#' @param ref_label Name of the reference checklist in the column Source.
#' @return A list of 5 data frames (Clade, Metric, Source), in the order
#'   Sper01, ITS2, 18S, Bota, reference.
family_richness_list <- function(df_H1, ref_label = "REF"){

  # Add the reference after the surveys
  tables <- c(
    df_H1$Result,
    list(tab_family_ref() %>% rename(Clade = family, Metric = total))
  )
  source_names <- c("Sper01", "ITS2", "18S", "Bota", ref_label)

  # Add the column Source
  lapply(seq_along(source_names), function(i) {
    df_temp <- tables[[i]]
    df_temp$Source <- source_names[i]
    return(df_temp)
  })
}

# ---- Growth form ----

#' Compute the share of each growth form in some families
#'
#' @param fam_vec Vector of families.
#' @return A data frame with family and one column per growth form (percentage
#'   of the species of the family in the reference checklist).
growth_form_family <- function(fam_vec){
  # Vector to df
  fam_df <- data.frame(family = fam_vec)

  # for Plot
  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ .x / rowSums(across(starts_with("stat."))) * 100)) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat.")) %>%
    inner_join(fam_df, by="family")
}

#' Estimate the share of each growth form in each survey
#'
#' The species of the reference checklist give the growth forms of each family.
#' They are weighted by the relative richness of the family in the survey.
#'
#' @param df_ab_rel Data frame with Clade (family), Metric (percentage of the
#'   survey) and Source (survey).
#' @return A data frame with Source and one column per growth form (percentage).
growth_form_source <- function(df_ab_rel){

  df_growth <- import_plants_nouragues() %>%
    rename(Clade = family)

  df_ab_rel %>%
    left_join( df_growth, by = "Clade") %>%
    mutate(across(starts_with("stat."), ~ .x * Metric)) %>%
    group_by(Source) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ .x / rowSums(across(starts_with("stat."))) * 100)) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat."))
}

#' Compute the share of each growth form in the reference checklist (Table 1)
#'
#' @return A data frame with type (growth form) and V1 (percentage), sorted by
#'   decreasing percentage.
tab_all_growth <- function(){

  import_plants_nouragues() %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ round(.x / rowSums(across(starts_with("stat."))) * 100, digits = 2))) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat.")) %>%
    t() %>%
    as.data.frame() %>%
    arrange(desc(V1)) %>%
    tibble::rownames_to_column("type")
}

#' Compute the share of each growth form in each family of the checklist
#'
#' @return A data frame with family and one column per growth form (percentage).
tab_family_growth <- function(){

  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(starts_with("stat."), ~ .x / rowSums(across(starts_with("stat."))) * 100)) %>%
    rename_with(~ sub("^stat\\.", "", .x), starts_with("stat."))
}

#' Compute the share of each family in the reference checklist
#'
#' @return A data frame with family and total (percentage of the species),
#'   sorted by decreasing percentage.
tab_family_ref <- function(){

  import_plants_nouragues() %>%
    group_by(family) %>%
    summarise(across(starts_with("stat."), ~ sum(.x, na.rm = TRUE)))%>%
    mutate(
      total = rowSums(across(starts_with("stat.")), na.rm = TRUE)
    ) %>%
    dplyr::select(family,total) %>%
    ungroup() %>%
    mutate(
      total = (total / sum(total)) * 100
    ) %>%
    arrange(desc(total))
}

# ===== 2. SPATIAL PATTERNS =====

# ---- Correlation ----

#' Correlate botanical and eDNA grids for each row of a table of parameters
#'
#' @param config_tab Data frame with the parameters of the two grids, see
#'   extract_configs(), and Metric_corr, see grid_correlation().
#' @return config_tab with Correlation, p.value and status. The rows without a
#'   correlation (too few observations, no variance) are removed.
make_cor_table <- function(config_tab){

  # Initialisation
  config_tab$Correlation <- vector("list", nrow(config_tab))

  # For each row in config tab
  config_tab_clean <- config_tab %>%
    mutate(
      Correlation = purrr::map(row_number(), ~ grid_correlation(config_tab[.x, ]))
    ) %>% # Correlation estimates
    tidyr::unnest_wider(Correlation) %>%  # r and p.value split
    dplyr::filter(!is.na(Correlation)) # r != NA

  return(config_tab_clean)
}

#' Correlate a botanical grid and an eDNA grid
#'
#' Metric_corr sets the method:
#'   "Count", "HYP", "EXP" : Spearman correlation between this model of the
#'                           botanical survey and the number of reads.
#'   "I.O"                 : Matthews correlation coefficient between the
#'                           presences of the two surveys.
#'   "OTU"                 : Spearman correlation between the metrics b.Metric
#'                           and g.Metric given in the row.
#'
#' @param row One row of the table of make_cor_table().
#' @return A list with Correlation and p.value, and status for Spearman.
grid_correlation <- function(row){

  # Split Metric_corr cases
  if(row$Metric_corr %in% c("Count","EXP","HYP")){
    row$b.Metric <- row$Metric_corr  # Bota
    row$g.Metric <- "Reads"          # Genet
  }
  else if(row$Metric_corr == "I.O"){
    row$b.Metric <- "Presence"   # Bota
    row$g.Metric <- "Presence"   # Genet
  }

  # Extract configs
  configs <- extract_configs(row)

  # Read Bota grid
  grid_bota <- grid_read(configs$bota)
  B <- grid_bota$df

  # Read Genet grid
  grid_genet <- grid_read(configs$genet)
  G <- grid_genet$df

  # Join
  BG <- inner_join(B, G, by = "sample_loc")

  # Select Final metric
  grid_corr <- switch(row$Metric_corr,
                      "Count" = get_spearman(B$Metric, G$Metric),
                      "HYP"   = get_spearman(B$Metric, G$Metric),
                      "EXP"   = get_spearman(B$Metric, G$Metric),
                      "OTU"   = get_spearman(B$Metric, G$Metric),

                      "I.O" = {
                        val <- mcc(BG$Metric.x, BG$Metric.y)
                        list(Correlation = val, p.value = NA)
                      },

                      stop(paste0("Unknown Metric : ", row$Metric_corr))
  )

  return(grid_corr)
}

#' Compute a Spearman correlation without failing on degenerate data
#'
#' @param x,y Numeric vectors.
#' @return A list with Correlation, p.value and status: "ok", "too_few_obs",
#'   "zero_variance" or "cor_test_failed" (Correlation is then NA).
get_spearman <- function(x, y){

  # Remove NA
  df <- na.omit(data.frame(x, y))

  # Not enough data
  if(nrow(df) < 3){
    return(list(
      Correlation = NA,
      p.value = NA,
      status = "too_few_obs"
    ))
  }

  # Zero variance
  if(sd(df$x) == 0 || sd(df$y) == 0){
    return(list(
      Correlation = NA,
      p.value = NA,
      status = "zero_variance"
    ))
  }

  # Safe test
  ct <- tryCatch(
    cor.test(df$x, df$y, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )

  if(is.null(ct)){
    return(list(
      Correlation = NA_real_,
      p.value = NA_real_,
      status = "cor_test_failed"
    ))
  }

  list(
    Correlation = unname(ct$estimate),
    p.value = ct$p.value,
    status = "ok"
  )
}

# ---- Other ----

#' Split a row of parameters into the configs of the two grids
#'
#' Columns prefixed with "b." go to the botanical grid, those prefixed with
#' "g." to the eDNA grid, and the others (Level, Name...) to both.
#'
#' @param row One row of a table of parameters.
#' @return A list with the configs bota and genet, see grid.R.
extract_configs <- function(row) {

  # Bota columns: prefixed with b. or shared
  b_cols <- grep("^(b\\.|(?!g\\.))", names(row), value = TRUE, perl = TRUE)
  bot_config <- as.list(row[, b_cols])
  bot_config$Type <- "bota"

  # Clean names: remove the prefix b.
  names(bot_config) <- sub("^b\\.", "", names(bot_config))

  # Genet columns: prefixed with g. or shared
  g_cols <- grep("^(g\\.|(?!b\\.))", names(row), value = TRUE, perl = TRUE)
  genet_config <- as.list(row[, g_cols])
  genet_config$Type <- "genet"

  # Clean names: remove the prefix g.
  names(genet_config) <- sub("^g\\.", "", names(genet_config))

  list(
    bota  = bot_config,
    genet = genet_config
  )
}

#' Compute the Matthews correlation coefficient between two presence vectors
#'
#' @param x,y Vectors of 0 and 1.
#' @return The coefficient, or NA if one margin of the confusion table is empty.
mcc <- function(x, y) {
  TP <- sum(x == 1 & y == 1, na.rm = TRUE)
  TN <- sum(x == 0 & y == 0, na.rm = TRUE)
  FP <- sum(x == 0 & y == 1, na.rm = TRUE)
  FN <- sum(x == 1 & y == 0, na.rm = TRUE)

  # Convert to numeric to avoid integer overflow
  TP <- as.numeric(TP)
  TN <- as.numeric(TN)
  FP <- as.numeric(FP)
  FN <- as.numeric(FN)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  if (is.na(denom) || denom == 0) return(NA_real_)

  (TP * TN - FP * FN) / denom
}

#' Summarise the morphological traits by taxon
#'
#' @param morpho Traits of the species, see import_plants_morpho().
#' @param level Name of the column with the taxonomic level.
#' @param fun Summary function, it must accept na.rm.
#' @return A data frame with the taxon, the summary of each numeric trait and n
#'   (number of species).
aggregate_morpho <- function(morpho, level = "Species", fun = median) {

  morpho %>%
    group_by(.data[[level]]) %>%
    summarise(
      across(
        where(is.numeric),
        ~ fun(.x, na.rm = TRUE)
      ),
      n = n(),
      .groups = "drop"
    )
}

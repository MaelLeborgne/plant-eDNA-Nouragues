
# ===== grid.R =====
# Builds the spatial grids of the plot: one value (Metric) for each of the
# 29 x 39 positions, spaced by 10 m, from a botanical or an eDNA survey.
#
# A grid is a list with two elements:
#   df     : data frame with x, y, sample_loc and Metric, one row per position.
#   config : list of the parameters used to build the grid.
#
# The parameters of a grid (config) are:
#   Type   : "bota" or "genet".
#   Level  : taxonomic level, "Family", "Genus" or "Species".
#   Name   : name of the taxon at this level, or "All" for every taxon.
#   Metric : value computed at each position.
#            bota  : "Count", "Presence", "OTU" (richness), "Shannon",
#                    "H1" (Hill number q=1), "HYP" and "EXP" (see below).
#            genet : "Reads", "Presence", "OTU" (richness), "Shannon", "H1".
#   Norm   : if TRUE, the values of the grid are divided by their sum.
#   Year   : bota only, inventory year (2012 or 2022).
#   Radius : bota only, distance in meters under which a stem is counted around
#            a position, or "pixel" for the 10 x 10 m square around it.
#   Alpha  : "HYP" and "EXP" only, exponent of the stem size.
#   Beta   : "EXP" only, decay rate with the distance.
#   Marker : genet only, "Sper01", "ITS2" or "18S".
#   Rar    : genet only, TRUE for the rarefied data.
#
# With d the distance between a stem and the position, and DBH = CircCorr / pi,
# the value at a position is the sum over the stems of:
#   HYP : DBH^Alpha / d
#   EXP : DBH^Alpha * exp(-Beta * d)
# The distance is bounded between 0.1 and 100 m.

# ===== Libraries =====
library(tidyverse) # Data manipulation
library(vegan)     # Diversity indices

# ===== INIT GRID =====

#' Create the positions of the plot with their soil sample
#'
#' @return A data frame with x, y, sample_loc (NA without a sample) and Metric
#'   set to zero.
ref_grid <- function(){

  # Create all positions
  df <- expand.grid(
    x = seq(10, 290, 10), # X position of samples
    y = seq(10, 390, 10)  # Y position of samples
  )

  # Import samples name reference
  samples_ref <- create_sample_ref()

  # Join samples name
  df <- left_join(df, samples_ref, by = c("x","y"))

  # Init Metric to zero
  df$Metric <- 0

  return(df)
}

#' Create an empty grid
#'
#' @param config List of the parameters of the grid, see the top of this file.
#' @return A grid with Metric set to zero.
init_grid <- function(config){

  # Empty list
  grid <- list()

  # Attach empty grid
  grid$df <- ref_grid()

  # Attach config
  grid$config <- config

  return(grid)
}

#' Divide the values of a grid by their sum
#'
#' @param df Data frame of a grid.
#' @return The data frame with the values of Metric summing to 1.
grid_normalize <- function(df){
  df$Metric <- df$Metric / sum(df$Metric)
  return(df)
}

# ==== TAXA FILTER ====

#' Keep the rows of a taxon
#'
#' @param Data Data frame with a column named Level.
#' @param Level Name of the column with the taxonomic level.
#' @param Name Name of the taxon, or "All" to keep every row.
#' @return The rows of Data for this taxon.
select_tax <- function(Data, Level, Name){

  if(Name == "All"){
    return(Data)
  }

  # Select the focus in Data
  Data[!is.na(Data[[Level]]) & Data[[Level]] == Name,]
}

# ===== BOTA =====

#' Count the stems of each taxon around each position
#'
#' @param G Grid with the config, see the top of this file.
#' @param B Botanical inventory, see import_bota().
#' @return A matrix with the positions in rows, the taxa of the level in columns
#'   and the number of stems in the cells.
bota_build_matrix <- function(G, B) {

  df <- G$df
  C  <- G$config
  R  <- C$Radius

  all_levels <- unique(B[[C$Level]])

  # Empty matrix
  mat <- matrix(0,
                nrow = nrow(df),
                ncol = length(all_levels))

  colnames(mat) <- all_levels

  for (pixel in 1:nrow(df)) {

    px <- df$x[pixel]
    py <- df$y[pixel]

    dij <- sqrt((px - B$Xfield)^2 +
                  (py - B$Yfield)^2)

    if (R == "pixel") {
      sel <- which(B$Xfield > px-5 &
                     B$Xfield < px+5 &
                     B$Yfield > py-5 &
                     B$Yfield < py+5)
    } else {
      sel <- which(dij < R & dij > 0)
    }

    if (length(sel) > 0) {
      counts <- table(B[[C$Level]][sel])
      mat[pixel, names(counts)] <- counts
    }
  }
  row.names(mat) <- df$sample_loc

  return(mat)
}

#' Find the stems around a position
#'
#' @param px,py Position in meters.
#' @param B Botanical inventory, see import_bota().
#' @param R Radius in meters, or "pixel" for the 10 x 10 m square.
#' @return A list with sel (indices of the selected stems) and dij (distance of
#'   every stem to the position, bounded between 0.1 and 100 m).
get_neighbors <- function(px, py, B, R){

  dij <- sqrt((px - B$Xfield)^2 +
                (py - B$Yfield)^2)

  dij[dij > 100] <- 100
  dij[dij < 0.1] <- 0.1

  if (R == "pixel") {
    sel <- which(B$Xfield > px-5 &
                   B$Xfield < px+5 &
                   B$Yfield > py-5 &
                   B$Yfield < py+5)
  } else {
    sel <- which(dij < R & dij > 0)
  }

  list(sel = sel, dij = dij)
}

#' Fill a grid from a botanical inventory
#'
#' @param grid Empty grid, see init_grid().
#' @return The grid with the Metric of its config at each position.
bota_grid_read <- function(grid) {

  # Configuration
  G <- grid
  C <- G$config
  M <- C$Metric
  R <- C$Radius

  # Importation
  bota <- import_bota(C)

  # Replace NA by Indet. (same as bota survey)
  for (l in c("Family","Genus","Species")) {
    bota[[l]][is.na(bota[[l]])] <- "Indet."
  }

  # Selection
  B <- select_tax(bota, C$Level, C$Name)

  # Build Matrix
  mat <- bota_build_matrix(G, B)

  # Fill Metric
  grid$df$Metric <- switch(M,
                           "Count"    = rowSums(mat),
                           "Presence" = as.integer(rowSums(mat) > 0),
                           "OTU"      = rowSums(mat > 0),
                           "Shannon"  = vegan::diversity(mat, index="shannon"),
                           "H1"       = exp(vegan::diversity(mat, index="shannon")),
                           "HYP"      = {
                             df <- G$df

                             sapply(1:nrow(df), function(pixel){

                               px <- df$x[pixel]
                               py <- df$y[pixel]

                               res <- get_neighbors(px, py, B, R=10000)
                               sel <- res$sel
                               dij <- res$dij

                               if (length(sel) == 0) return(0)

                               sum(((B$CircCorr[sel]/pi)^C$Alpha) / dij[sel],
                                   na.rm = TRUE)
                             })
                           },
                           "EXP" = {
                             df <- G$df

                             sapply(1:nrow(df), function(pixel){

                               px <- df$x[pixel]
                               py <- df$y[pixel]

                               res <- get_neighbors(px, py, B, R)
                               sel <- res$sel
                               dij <- res$dij

                               if (length(sel) == 0) return(0)

                               sum(((B$CircCorr[sel]/pi)^C$Alpha) *
                                     exp(-C$Beta * dij[sel]),
                                   na.rm = TRUE)
                             })
                           },
                           stop("Unknown Metric: ", M)
  )

  # Normalization
  if(isTRUE(C$Norm)){
    grid$df <- grid$df %>% grid_normalize()
  }

  # Return
  return(grid)
}

# ==== GENET READ ====

#' Fill a grid from the eDNA data of a marker
#'
#' Positions without a soil sample keep the value zero.
#'
#' @param grid Empty grid, see init_grid().
#' @return The grid with the Metric of its config at each position.
genet_grid_read <- function(grid){

  # Configuration
  G <- grid
  C <- G$config
  M <- C$Metric

  # Import genet
  genet <- import_genet(C)

  # Replace NA by Indet. (same as bota survey)
  for (l in c("family","genus","species")) {
    genet$motus[[l]][is.na(genet$motus[[l]])] <- "Indet."
  }

  # Select OTU
  otu_list <- select_tax(genet$motus, tolower(C$Level), C$Name)

  # get otu_id
  otu_list <- otu_list %>% filter(!is.na(id)) %>% pull(id)

  # Replace : and - by . # (18S anomaly)
  otu_list <- gsub(":", ".", otu_list)
  otu_list <- gsub("-", ".", otu_list)

  # Select Reads from OTU_id
  reads_filtered <- genet$reads[, colnames(genet$reads) %in% otu_list, drop = FALSE]

  # Init
  reads_tax <- reads_filtered

  # Only for Shannon and H1 above the species level:
  # group the OTU of the same taxon
  if(M %in% c("Shannon","H1") & C$Level != "Species"){

    # Taxon of each OTU, aligned with the columns
    tax_names <- genet$motus[[tolower(C$Level)]][
      match(colnames(reads_filtered), genet$motus$id)
    ]
    # Manage NA
    valid <- !is.na(tax_names)

    # Sum the reads by taxon
    reads_tax <- t(rowsum(
      t(reads_filtered[, valid]),
      group = tax_names[valid]
    ))
  }

  # Create dataframe for calcul
  metric_df <- data.frame(sample_loc = rownames(reads_tax))

  # Compute the metric
  metric_df$Metric <- switch(M,
                    "Reads"         = rowSums(reads_tax),
                    "OTU"           = rowSums(reads_tax > 0),
                    "Presence"      = as.integer(rowSums(reads_tax > 0) > 0),
                    "Shannon"       = diversity(reads_tax, index = "shannon"),
                    "H1"            = exp(diversity(reads_tax, index = "shannon")),
                    stop("Unknown Metric: ", M)
  )

  # Join by sample_loc
  grid$df <- left_join(grid$df, metric_df, by="sample_loc")

  # Merge Metric cols
  grid$df <- grid$df %>%
    mutate(Metric = coalesce(Metric.y, Metric.x)) %>%
    dplyr::select(-Metric.x, -Metric.y)

  # Normalization
  if(isTRUE(C$Norm)){
    grid$df <- grid$df %>% grid_normalize()
  }

  return(grid)
}

# ==== GRID READ ====

#' Build a grid from its parameters
#'
#' @param config List of the parameters of the grid, see the top of this file.
#' @return The filled grid.
grid_read <- function(config){

  # Initialize Grid
  grid <- init_grid(config)

  # Fill Grid
  if(config$Type == "bota"){
    grid <- bota_grid_read(grid)
  }
  else if (config$Type == "genet"){
    grid <- genet_grid_read(grid)
  }
  else{
    stop(paste0("Unknown config$Type : ", config$Type))
  }

  return(grid)
}

#' Build one grid per row of a table of parameters
#'
#' @param config_tab Data frame with one column per parameter of the grid.
#' @return A list of filled grids, in the order of the rows.
build_grids <- function(config_tab){

  # Iterate on rows
  grid_list <- lapply(seq_len(nrow(config_tab)), function(i){
    # get parameters
    config <- as.list(config_tab[i, ])
    # Import grid
    grid_read(config)
  })

  return(grid_list)
}

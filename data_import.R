
# ===== data_import.R =====
# Paths of the project and functions reading the raw data files.
# Every other module gets its data through the import_*() functions below.

# ===== Libraries =====
library(tidyverse) # Data manipulation
library(readxl)    # Read Excel files

# ===== Paths =====

# Root of the project. Files are read relatively to the working directory,
# which is the project folder when app.Rmd is knitted or opened in RStudio.
PROJECT_PATH <- getwd()

# Raw data files
DATA_PATH <- file.path(PROJECT_PATH, "data")

# Intermediate results that are long to compute
CACHE_PATH <- file.path(PROJECT_PATH, "ADNe_Nourague_2026_files")

# ===== Cache =====

#' Load an intermediate result from the cache, or compute and save it
#'
#' @param name Name of the result, the file is CACHE_PATH/.<name>.RData.
#' @param compute Function without argument that computes the result.
#' @param recompute If TRUE, compute the result even if it is in the cache.
#' @return The result.
cached <- function(name, compute, recompute = FALSE){

  file <- file.path(CACHE_PATH, paste0(".", name, ".RData"))

  # Read the saved result
  if(!recompute && file.exists(file)){
    env <- new.env()
    load(file, envir = env)
    return(env[[name]])
  }

  # Compute and save, the object has the same name as the result
  env <- new.env()
  env[[name]] <- compute()
  save(list = name, envir = env, file = file)

  return(env[[name]])
}

# ===== SAMPLES POSITION (L.Zinger) =====

#' Read the position of the soil samples in the plot
#'
#' @return A data frame with one row per sampling location:
#'   sample_loc (location name), x and y (position in meters).
create_sample_ref <- function(){

  # Read TXT
  data_samples <- read.table(file.path(DATA_PATH, "Zinger_plants_nouragues",
                                       "Zinger_et_al_metabar_PP1_samples.txt"),
                             header = TRUE,
                             comment.char = "#",
                             fill = TRUE,
                             blank.lines.skip = TRUE)

  # Extracts samples name and position
  sample_ref <- data_samples %>%
    dplyr::select(pcr,x,y) %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    distinct(sample_loc, x, y)

  return(sample_ref)
}

# ==== BOTA (G.Sellan & A.Mirabelle) ====

#' Read a botanical inventory of the plot
#'
#' @param config List with the element Year: inventory year, 2012 or 2022.
#' @return A data frame with one row per stem. Main columns: Family, Genus,
#'   Species (binomial name), Xfield and Yfield (position in meters) and
#'   CircCorr (circumference of the stem).
import_bota <- function(config){

  # Select Year
  if(config$Year == "2012"){
    file <- "Nour_2012.csv"
  }
  else if(config$Year == "2022"){
    file <- "Nour_2022.csv"
  }
  else{
    stop("Unknown config$Year : ", config$Year)
  }

  read.csv(file.path(DATA_PATH, "Bota", file), sep = ",", dec = ".") %>%
    mutate(Species = paste(Genus, Species, sep=" "))
}

# ==== GENET (I.Calderon) ====

#' Read the clustered eDNA data of a marker, before rarefaction
#'
#' @param config List with the element Marker: "Sper01", "ITS2" or "18S".
#' @return A metabaR list, see import_genet().
select_clust <- function(config){

  file <- switch(config$Marker,
    "Sper01" = "Sper01_mblist_annot_metabar_clean_clust97.rds",
    "18S"    = "18S_mblist_annot_metabar_clean.rds",
    "ITS2"   = "ITS2_mblist_annot_metabar_clean_clust97.rds",
    stop("Unknown Marker: ", config$Marker)
  )

  readRDS(file.path(DATA_PATH, "Irene_genet", file))
}

#' Read the clustered and rarefied eDNA data of a marker
#'
#' @param config List with the element Marker: "Sper01", "ITS2" or "18S".
#' @return A metabaR list, see import_genet().
select_rar <- function(config){

  file <- switch(config$Marker,
    "Sper01" = "Sper01_mblist_annot_metabar_clean_clust97_RAR.rds",
    "18S"    = "18S_mblist_annot_metabar_clean_clust97_RAR.rds",
    "ITS2"   = "ITS2_mblist_annot_metabar_clean_clust97_RAR.rds",
    stop("Unknown Marker: ", config$Marker)
  )

  readRDS(file.path(DATA_PATH, "Irene_genet", file))
}

#' Read the eDNA data of a marker
#'
#' @param config List with the elements Marker ("Sper01", "ITS2" or "18S") and
#'   Rar (TRUE for rarefied data, all the figures of the report use them).
#' @return A metabaR list. Main elements: reads (matrix of read counts, samples
#'   in rows and MOTUs in columns) and motus (data frame with the id and the
#'   taxonomic assignment of each MOTU: family, genus, species).
import_genet <- function(config){

  if(config$Rar){
    return(select_rar(config))
  }

  select_clust(config)
}

# ==== ABIOTIC (L.Zinger) ====

#' Read the environmental variables of the sampling locations
#'
#' @return A data frame with one row per location: x and y (position in meters)
#'   and the edaphic and light variables (wetness, Al, K, Ca, Na, Fe, C, N,
#'   light...).
import_abiotic <- function(){

  # Read TXT
  data <- read.table(file.path(DATA_PATH, "Zinger_plants_nouragues",
                               "Zinger_et_al_metabar_PP1_samples.txt"),
                     header = TRUE,
                     comment.char = "#")

  # Split samples name
  abiotic <- data %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    filter(!control) %>% # Remove controls
    dplyr::select(-c(tag_combo, control,sample, sample_loc, sample_rep, chempred)) %>% # select
    dplyr::distinct(x, y, .keep_all = TRUE) # Remove replicats

  return(abiotic)
}

# ==== PLANTS MORPHO (Vleminckx et al. 2021, provided by G.Sellan) ====

#' Read the morphological traits of the species
#'
#' @return A data frame with the taxonomy (Family, Genus, Species) and one
#'   numeric column per trait, sheet "App.S6-ok" of Appendix_S2-6___S8-13.xlsx.
import_plants_morpho <- function(){

  # Read XLS
  data <- read_excel(file.path(DATA_PATH, "Giacomo_plant_morpho",
                               "Appendix_S2-6___S8-13.xlsx"),
                     sheet = "App.S6-ok", skip = 3) %>%
    mutate(Species = gsub("_", " ", Species))

  return(data)
}

# ==== PLANTS Nouragues (Bongers et al. 2001, provided by L.Zinger) ====

#' Read the reference checklist of the plants of the Nouragues high-forest
#'
#' @return A data frame with one row per species: family, genus, species and
#'   the growth forms as 0/1 columns prefixed with "stat." (stat.tree...).
import_plants_nouragues <- function(){

  # Read XLS
  data <- read_excel(file.path(DATA_PATH, "Zinger_plants_nouragues",
                               "METABAR_plantes", "metadata",
                               "Plants_Nouragues2001.xls")) %>%
    dplyr::filter(habit.highforest == 1)

  return(data)
}

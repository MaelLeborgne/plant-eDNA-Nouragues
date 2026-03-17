
# ===== data_import.R =====
# Regroup importation functions by Subject

# ===== Libraries =====
library(tidyverse) # Basic
library(dplyr)     # Basic
library(readxl)    # Excel

# ===== Source =====
PROJECT_PATH = "/home/jean/Bureau/ADNe_Nourague_2026/"


# ===== SAMPLES POSITION (L.Zinger) =====

create_sample_ref <- function(){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Zinger_plants_nouragues/"))
  
  # Read TXT
  data_samples <- read.table("Zinger_et_al_metabar_PP1_samples.txt",
                             header = TRUE,
                             comment.char = "#")
  
  # Extracts samples name and position
  sample_ref <- data_samples %>%
    select(pcr,x,y) %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    distinct(sample_loc, x, y)
  
  return(sample_ref)
}

# ==== BOTA (G.Sellan & A.Mirabelle) ====

import_bota <- function(config){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Bota/"))
  
  # Select Year
  if(config$Year == "2012"){
    bota <- read.csv("Nour_2012.csv",sep = ",",dec = ".") %>%
      mutate(Species = paste(Genus, Species, sep=" "))
  }
  else if(config$Year == "2022"){
    bota <- read.csv("Nour_2022.csv",sep = ",",dec = ".") %>%
      mutate(Species = paste(Genus, Species, sep=" "))
  }
  else{
    stop("Unknown Year: ", config$Year)
  }

  return(bota)
}

# ==== GENET (I.Calderon) ====


# Data clustered
select_clust <- function(config){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Irene_genet/"))
  
  # Select Marker
  data <- switch (config$Marker,
    "GH" = readRDS("GH_mblist_annot_metabar_clean_clust97.rds"),
    "18S"= readRDS("18s_mblist_annot_metabar_clean.rds"),
    "ITS"= readRDS("GH_mblist_annot_metabar_clean_clust97.rds"),
    stop("Unknown Marker: ", config$Marker)
  )
  
  return(data)
}

# Data rarefied
select_rar <- function(config){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Irene_genet/"))
  
  # Select Marker
  data <- switch (config$Marker,
                  "GH" = readRDS("GH_mblist_annot_metabar_clean_clust97_RAR.rds"),
                  "18S"= readRDS("18S_mblist_annot_metabar_clean_clust97_RAR.rds"),
                  "ITS"= readRDS("ITS_mblist_annot_metabar_clean_clust97_RAR.rds"),
                  stop("Unknown Marker: ", config$Marker)
  )
  
  return(data)
}

# Apply option of Rarification
import_genet <- function(config){
  
  # Select rarefied data
  if(config$Rar){ 
    data <- select_rar(config)
  }
  else{ # Or not
    data <- select_clust(config)
  }
  
  return(data)
}

# ==== ABIOTIC (L.Zinger) ====

import_abiotic <- function(){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"Zinger_plants_nouragues/"))
  
  # Read TXT
  data <- read.table("Zinger_et_al_metabar_PP1_samples.txt",
                             header = TRUE,
                             comment.char = "#")
  
  # Split samples name
  abiotic <- data %>%
    separate(pcr, into = c("sample_loc", "sample_rep"), sep = "_") %>%
    filter(!control) # Remove controls
  
  return(abiotic)
}

# ==== PLANTS MORPHO (Giacomo) ====

import_plants_morpho <- function(){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Giacomo_plant_morpho/"))
  
  # Read XLS
  data <- read_excel("Appendix_S2-6___S8-13.xlsx", sheet = "App.S6-ok", skip = 3)
  
  return(data)
}

# ==== PLANTS Nouragues (L.Zinger) ====

import_plants_nouragues <- function(){
  
  # Set directory
  setwd(paste0(PROJECT_PATH,"data/Zinger_plants_nouragues/METABAR_plantes/metadata/"))
  
  # Read XLS
  data <- read_excel("Plants_Nouragues2001.xls")
  
  return(data)
}


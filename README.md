# Fine-scale distribution of plant eDNA in a tropical rainforest

R code behind the figures of the master's thesis *Fine-scale distribution of plant eDNA in a tropical rainforest* (Maël Leborgne, Université des Antilles, UMR EcoFog, 2026), supervised by Irène Calderón-Sanou, Lucie Zinger and Jérémy Raynaud.

The study compares two views of the plants of the Petit Plateau plot (Nouragues reserve, French Guiana): a botanical inventory of the stems, and the plant DNA found in soil samples taken every 10 m, with three genetic markers (Sper01, ITS2, 18S). It looks at the taxa detected by each survey, how well their spatial patterns match, and which environmental conditions drive their diversity.

## Contents

| File | Content |
|---|---|
| `app.Rmd` | **Start here.** Generates each figure of the report, in the order of the report, with its explanation |
| `data_import.R` | Paths of the project and functions reading the raw data |
| `grid.R` | Spatial grids of the plot, with the description of all their parameters |
| `table.R` | Taxa lists, relative richness, growth forms and correlations |
| `plots.R` | Functions building the figures |

## Data

**The data are not in this repository.** They belong to the people and projects that produced them. Please cite them and respect their conditions of use.

| Data | Files | Producers | Access |
|---|---|---|---|
| Soil samples and environmental variables | `Zinger_et_al_metabar_PP1_samples.txt` | Zinger et al. (2019) | Public, [Dryad](https://doi.org/10.5061/dryad.dc28m) |
| Reference checklist of the plants | `Plants_Nouragues2001.xls` | From Bongers et al. (2001), provided by L. Zinger | Contact L. Zinger |
| Morphological traits | `Appendix_S2-6___S8-13.xlsx` | Vleminckx et al. (2021) | Public, [Dryad](https://doi.org/10.5061/dryad.c59zw3r76) (CC0) |
| Botanical inventories 2012 and 2022 | `Nour_2012.csv`, `Nour_2022.csv` | Guyafor network, access granted by A. Mirabelle | Contact the data owners |
| Cleaned eDNA data | `*_mblist_annot_metabar_clean*.rds` | L. Zinger, cleaned by A. Maniouloux and I. Calderón-Sanou | Contact the data owners |

The program expects a folder `Data-Nouragues-2026` next to the folder of the program:

```
parent_folder/
├── plant-eDNA-Nouragues/              (this repository)
└── Data-Nouragues-2026/
    ├── Bota/
    │   ├── Nour_2012.csv
    │   └── Nour_2022.csv
    ├── Giacomo_plant_morpho/
    │   └── Appendix_S2-6___S8-13.xlsx
    ├── Irene_genet/
    │   ├── 18S_mblist_annot_metabar_clean.rds
    │   ├── 18S_mblist_annot_metabar_clean_clust97_RAR.rds
    │   ├── ITS2_mblist_annot_metabar_clean_clust97.rds
    │   ├── ITS2_mblist_annot_metabar_clean_clust97_RAR.rds
    │   ├── Sper01_mblist_annot_metabar_clean_clust97.rds
    │   └── Sper01_mblist_annot_metabar_clean_clust97_RAR.rds
    └── Zinger_plants_nouragues/
        ├── Zinger_et_al_metabar_PP1_samples.txt
        └── Plants_Nouragues2001.xls
```

## Usage

The code was last run with R 4.6.1. Install the packages:

```r
install.packages(c(
  "tidyverse", "readxl", "vegan", "ggtext", "ggpubr", "viridis", "patchwork",
  "akima", "emmeans", "multcomp", "multcompView", "randomForest",
  "ggVennDiagram", "UpSetR", "car", "FSA", "rcompanion", "scales", "rmarkdown"
))
```

Then open `ADNe_Nourague_2026.Rproj` in RStudio and knit `app.Rmd`, or from the folder of the program:

```r
rmarkdown::render("app.Rmd")
```

The first run computes all the results from the data, which took about 15 minutes on a laptop. They are then saved in `ADNe_Nourague_2026_files/` and read from there. To compute them again, set `RECOMPUTE <- TRUE` at the top of `app.Rmd`.

## Reproducibility

Recomputing from the data gives the same results as the report, except for the random forests of Figure 11, which vary slightly because they are random (their median importances differed by less than 0.003 in our test).

## References

- Bongers, F., Charles-Dominique, P., Forget, P.-M., & Théry, M. (Eds). (2001). *Nouragues*. Springer Netherlands. https://doi.org/10.1007/978-94-015-9821-7
- Vleminckx, J., Fortunel, C., Valverde-Barrantes, O., et al. (2021). Resolving whole-plant economics from leaf, stem and root traits of 1467 Amazonian tree species. *Oikos*, 130(7), 1193–1208. https://doi.org/10.1111/oik.08284
- Zinger, L., Taberlet, P., Schimann, H., et al. (2019). Body size determines soil community assembly in a tropical forest. *Molecular Ecology*, 28(3), 528–543. https://doi.org/10.1111/mec.14919

## Acknowledgments

This work relies on the people who produced the botanical and genetic data, from plant identification in the field to the sequencing and processing of the genetic data. See the acknowledgments of the thesis.

## License

No license has been chosen yet: all rights are reserved until one is added. The data are not covered by any license of this repository.

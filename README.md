# KenyaAtlasComparison

This repository contains the code and data used for the following publication:

> Nussbaumer, R., ....

The main output of this code can be visualized on [kenyabirdtrends.co.ke](https://kenyabirdtrends.co.ke/)

[![image](https://github.com/user-attachments/assets/7ebb2728-29a3-4171-ad4a-60d4aa79ac33)](https://kenyabirdtrends.co.ke/)
_Explore The change in distributation of any species in the species view and searching for you species of interest_


## Code structure

### Core analysis

You can reproduce the analysis by running the script in alphabetical order:

1. `A_import_old_atlas.m`: Process historical data (`data/oldatlas/A Bird Atlas of Kenya_v5.xlsx`) based on the species list `data/species_base_list.xlsx` into `data/oldatlas.mat` as well as the grid data (read from `data/oldatlas/grid.geojson`)
2. `B_import_KBM.m`: Process KBM data (downloaded from the ABAP API in `data/kbm/geojson/`) based on the species list `data/kbm/sp_kbm.xlsx` into `data/kbmatlas.mat`
3. `C_import_ebird.m`: Process eBird data (EBD file downloaded in `data/eBird/ebd_KE_relOct-2023/`) based on the species list `data/eBird/sp_ebird.xlsx` into `data/ebirdatlas.mat` 
4. `D_correction.m`: Compute the confidence index into the corrected grid `data/grid_corr.mat`

### Export

- `E_export_figure.m`: Generate figures for the paper
- `E_export_website.m`: Generate dataset in `export/website/` for [kenyabirdtrends.co.ke](https://kenyabirdtrends.co.ke/)

### Post-analysis for paper

`F_analysis.m`: More figures, some used in the paper, other not.
`F_phylo.m`: Phyologenetic data

### Other analysis and comparison not used in paper

- `G_burns.m`: Comparison with [Burns et al. (2021)](https://doi.org/10.1002/ece3.8282)
- `G_CSR8.m`: Comparison with the [CSR 8 population trends from AEWA](https://iwc.wetlands.org/index.php/aewatrends8)
- `G_Ngulia.m`: Comparison with [Pearson et al. (2017)](https://www.ajol.info/index.php/scopus/article/view/149917)
- `G_sabap.m`: Comparison with [Underhill & Brooks (2014)](https://journals.uct.ac.za/index.php/BO/article/view/235)
- `G_scagenger.m`: 
- `G_serengeti.m`: Comparison with [Henao-Diaz & Sinclair (2019)](https://doi.org/10.1002/ecy.2919)
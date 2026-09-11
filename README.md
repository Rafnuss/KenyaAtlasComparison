[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13710126.svg)](https://doi.org/10.5281/zenodo.13710126)

# KenyaAtlasComparison

This repository contains the code and data used for the following publication:

> Nussbaumer, R., Nussbaumer, A., Guchu, S., Hatfield, R.S., M. Kanga, E., Kung'u, G.N., Kuria, A., Miller, E., Ndang'ang'a, P.K., Njoroge, P., Ogada, D., Shema, S. and Jackson, C. (2024), Historical Bird Atlas and Contemporary Citizen Science Data Reveal Long-Term Changes in Geographic Range of Kenyan Birds. _Diversity and Distributions_ e13935. https://doi.org/10.1111/ddi.13935

The main output of this code can be visualized on [kenyabirdtrends.co.ke](https://kenyabirdtrends.co.ke/)

[![image](https://github.com/user-attachments/assets/7ebb2728-29a3-4171-ad4a-60d4aa79ac33)](https://kenyabirdtrends.co.ke/)
_Explore The change in distributation of any species in the species view and searching for you species of interest_


## Code structure

### Core analysis

You can reproduce the analysis by running the script in alphabetical order:

1. `A_import_old_atlas.m`: Process historical data (`data/oldatlas/A Bird Atlas of Kenya_v5.xlsx`) based on the species list `data/species_base_list.csv` into `data/oldatlas.mat` as well as the grid data (read from `data/oldatlas/grid.geojson`)
2. `B_import_KBM.m`: Process KBM data (downloaded from the ABAP API in `data/kbm/geojson/`) based on the species list `data/kbm/sp_kbm.xlsx` into `data/kbmatlas.mat`
3. `C_import_ebird.m`: Process eBird data (EBD file downloaded in `data/eBird/ebd_KE_relOct-2023/`) based on the species list `data/eBird/sp_ebird.xlsx` into `data/ebirdatlas.mat` 
4. `D_correction.m`: Compute the confidence index into the corrected grid `data/grid_corr.mat`

### Taxonomy

`data/species_base_list.csv` is the master species list. Each row is one
**atlas concept**, keyed by the immutable `SEQ` used since the 1989 atlas —
and a concept is not always one species: 35 of them are lumps that the atlas
recorded at a coarser resolution than today's checklists.

`data/taxonomy/update_taxonomy.py` resolves every concept against
[AviList](https://www.avilist.org) v2025b and the eBird/Clements taxonomy,
adding an `avibase_id` plus current names, family, IUCN status and eBird codes
to each row. `data/eBird/add_avibase_id.py` does the equivalent for
`sp_ebird.xlsx`, so eBird observations join on a stable concept id rather than
on a scientific name that changes with every taxonomic revision. Both scripts
are idempotent. See [`data/taxonomy/SOURCES.md`](data/taxonomy/SOURCES.md) for
the reference versions, the column contract and the conventions.

### Export

- `E_export_figure.m`: Generate figures for the paper
- `E_export_website.m`: Generate dataset in `export/website/` for [kenyabirdtrends.co.ke](https://kenyabirdtrends.co.ke/)

`export/website/sp_base.json` is a contract with
[Rafnuss/KenyaBirdTrends](https://github.com/Rafnuss/KenyaBirdTrends): pull it
in there with `npm run sync:data`, which validates the columns the site reads
(`tests/dataSchema.test.js` asserts the same set). Each species carries both
namings — the 1989 atlas one and the current AviList one — so the site can
offer them as alternatives; see the taxonomy section above.

### Post-analysis for paper

`F_analysis.m`: More figures, some used in the paper, other not.
`F_phylo.m`: Phyologenetic data. Needs the Bioinformatics Toolbox for
`phytreeread`; without it the script cannot run at all.

### Other analysis and comparison not used in paper

- `G_burns.m`: Comparison with [Burns et al. (2021)](https://doi.org/10.1002/ece3.8282)
- `G_CSR8.m`: Comparison with the [CSR 8 population trends from AEWA](https://iwc.wetlands.org/index.php/aewatrends8)
- `G_Ngulia.m`: Comparison with [Pearson et al. (2017)](https://www.ajol.info/index.php/scopus/article/view/149917)
- `G_sabap.m`: Comparison with [Underhill & Brooks (2014)](https://journals.uct.ac.za/index.php/BO/article/view/235)
- `G_scavenger.m`: Vulture trends against 1980 vs 2000 density estimates per season (hard-coded in the script, source not recorded)
- `G_serengeti.m`: Comparison with [Henao-Diaz & Sinclair (2019)](https://doi.org/10.1002/ecy.2919)

The SABAP1-2, Burns et al. and AEWA files were keyed to the atlas by hand,
before the 2026-09 taxonomy update retired some `SEQ`s through `merged_SEQ`.
Any script joining such a file on `SEQ` must first call
`functions/fold_SEQ.m`, which follows the taxonomic folds and zeroes the
rest (so the usual `t = t(t.SEQ>0,:)` filter drops them) — otherwise those
rows silently vanish from the comparison. `G_scavenger.m` and
`G_Ngulia.m` key on hard-coded `SEQ`s/names that no fold touches, and
`G_serengeti.m` does not use the atlas list at all.

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

Run `AA_wrapper_all.m` from this repository's root. It prepares taxonomy metadata and runs A–E in order:

1. `A_import_old_atlas.m`: Read the sibling digitization project's `export/comparison/atlas_records.csv`, select 1970–1984 records, and apply the explicit historical-to-comparison mapping.
2. `B_import_KBM.m`: Join cached KBM observations by source `ref` through `data/kbm/species_mapping.csv`.
3. `C_import_ebird.m`: Join archived eBird observations through `data/eBird/species_mapping.csv`.
4. `D_correction.m`: Refit the existing effort correction and export species statistics.
5. `E_export_website.m`: Export Avibase-keyed maps, original names, AviList names and associated source links.

### Taxonomy

Four readable CSVs contain the curated decisions: `data/species_base_list.csv`, `data/taxonomy/atlas_mapping.csv`, and the two source mapping files above. Each comparison taxon has one unique Avibase ID and explicit AviList membership. Historical concepts remain owned by `Digitization-of-A-Bird-Atlas-of-Kenya`; comparison merges and exclusions are recorded separately rather than changing that historical taxonomy.

`python3 data/taxonomy/update_taxonomy.py` derives display metadata into `data/taxonomy/comparison_taxa.csv` from those memberships and the digitization project's pinned AviList 2026 export. eBird observations remain a data source; eBird and Birds of Kenya 2019 naming systems are retired. Original and AviList names travel together to the website. See [taxonomy sources and decisions](data/taxonomy/SOURCES.md) for the table contracts, qualifications and rebuild steps.

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

The SABAP1-2, Burns et al. and AEWA files retain legacy numeric SEQs. `functions/fold_SEQ.m` reads the explicit `atlas_mapping.csv` decisions before these scripts join them to current comparison taxa. Identification merges are not applied to external-study trends. `G_scavenger.m` and `G_Ngulia.m` use fixed legacy SEQs; `G_serengeti.m` does not use the atlas list.

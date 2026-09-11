# Taxonomic reference sources

Reference files used to build the `SEQ` <-> AviList crosswalk that lives directly
in `data/species_base_list.csv` (see below). Re-download with
`bash data/taxonomy/fetch_sources.sh`. **Pin these versions**: the crosswalk
columns are only meaningful against the versions listed here.

| File | Version | Retrieved | Rows | Key columns |
|---|---|---|---|---|
| `AviList-v2025b-10Jun2026-extended.xlsx` | v2025b (10 Jun 2026) | 2026-09-10 | 33,686 taxa / 11,131 species | `Scientific_name`, `English_name_AviList`, `AvibaseID`, `Species_code_Cornell_Lab`, `Order`, `Family_English_name`, `IUCN_Red_List_Category` |
| `AviList_v2025_metadata_11Jun.pdf` | v2025 | 2026-09-10 | — | column definitions |
| `eBird_taxonomy_v2025-4.csv` | eBird taxonomy v2025.4 | 2026-09-10 | 17,891 | `SPECIES_CODE`, `SCI_NAME`, `PRIMARY_COM_NAME`, `CATEGORY`, `REPORT_AS`, **`TAXON_CONCEPT_ID`** |
| `avilist_v2025b_species.csv` | derived | 2026-09-10 | 11,131 species | subset of the extended xlsx, species rank only - regenerate from the xlsx if columns are ever added |
| `ebd_KE_relOct-2023_records_per_concept.tsv` | derived from `ebd_KE_relOct-2023.zip` | 2026-09-10 | 1,573 | `TAXON_CONCEPT_ID`, `n_records` - Kenyan EBD record counts per taxon, used only to document historical-lump comments (see `update_taxonomy.py`) |

## Why two taxonomies

- **AviList v2025b** is the target taxonomy: its extended sheet carries both
  `AvibaseID` and `Species_code_Cornell_Lab` in the same row, bridging
  `species_base_list.csv` to eBird in one step.
- **eBird taxonomy** is needed because the raw EBD reports `slash`, `spuh`,
  `issf` and `hybrid` taxa that AviList (species + subspecies only) does not
  contain, and `eBird_taxonomy_v2025-4.csv` specifically is the only source
  that carries an Avibase id (`TAXON_CONCEPT_ID`) for those broader taxa -
  essential for resolving lumped atlas concepts. The API variant
  (`api.ebird.org/v2/ref/taxonomy/ebird`) lacks that column and was dropped.

`AviList-v2025b-10Jun2026-short.xlsx` was removed - it is a column subset of
the extended file with nothing the extended file doesn't have.

## Manual step

Cloudflare blocks scripted access to `birds.cornell.edu`, so
`eBird_taxonomy_v2025-4.csv` must be downloaded by hand - see
`fetch_sources.sh` for the exact steps. The eBird REST API
(`https://ebird.org/api/keygen` for a free key) is scriptable but does not
carry `TAXON_CONCEPT_ID`, so it cannot substitute for this file.

## Citations

> AviList Core Team. 2026. AviList: The Global Avian Checklist, v2025b.
> <https://doi.org/10.2173/avilist.v2025b>. CC BY 4.0.

> Clements, J. F., P. C. Rasmussen, T. S. Schulenberg, et al. 2025. The eBird/Clements
> Checklist of Birds of the World: v2025. <https://www.birds.cornell.edu/clementschecklist/>

## Next AviList release

AviList updates annually (~June); the next release is expected autumn 2026,
and the following eBird/Clements update is intended to be *fully* aligned with
it. v2025b and eBird/Clements v2025 are currently 99.993% aligned (a few dozen
species differ) - see `update_taxonomy.py`'s `avibase_id_todo.csv` output for
any that surface as unresolved in this dataset.

---

## The crosswalk: one file, not two

Step 3 (linking each `SEQ` to its current AviList species) is **not** a
separate join table. Every atlas concept is a lump of 1+ modern species, and
`data/eBird/sp_ebird.xlsx` already records that membership (it maps every
eBird taxon, including slash/spuh, onto a `SEQ`) - so no new many-to-many file
is needed. `update_taxonomy.py` instead adds derived columns directly onto
`data/species_base_list.csv`:

| Column | Content |
|---|---|
| `avilist_common_name` | display name: the eBird slash's own name when one exactly covers the concept, else a generated slash label, else the single AviList English name |
| `avilist_scientific_name` | single displayable binomial: the species' own when there is one, else the shared genus + `sp.` (`Ficedula sp.`), else the members joined with `/` |
| `avilist_members` | every member binomial, `\|`-joined - the precise content of the concept |
| `avilist_sort` | AviList's own linear sequence, for ordering a species list "by taxonomy"; a lump takes its earliest member's position |
| `ebird_code` | member eBird species code(s), `\|`-joined |
| `family` | AviList English family name (from the member(s)), falling back to eBird's own family when AviList has no species-rank row - `SEQ 198` (African Swamphen, which AviList ranks as a subspecies of Purple Swamphen while eBird splits it) and `SEQ 556` (resolvable only to a slash) would otherwise be blank |
| `order` | AviList order (from the member(s)), same fallback |
| `iucn` | IUCN Red List code. For a lump this is the **most severe** status among its members, the usual conservation-reporting convention: `SEQ 753` (Bar-throated/Taita Apalis) contains Critically Endangered *Apalis fuscigularis*, and reporting the concept as blank would drop it out of every by-threat-category analysis. 5 concepts are affected. |
| `birdlife_url` | AviList's BirdLife DataZone factsheet, only when `n_avilist_species == 1` - a lump has no single factsheet to point at |
| `n_avilist_species` | count of distinct modern species making up the concept (1 for the great majority; up to 5 for `SEQ 786`) |

MATLAB reads a lump's pipe-delimited columns with `strsplit(x, "|")`.

Column names are snake_case throughout, and the `avilist_*` triplet
(`_common_name`, `_scientific_name`, `_sort`) deliberately mirrors the
historical `common_name` / `scientific_name` / `SEQ` triplet: the website
offers the two as alternative namings of the same row, so they are addressed
the same way. See `TAXONOMY_FIELDS` in Rafnuss/KenyaBirdTrends `src/store.js`.

## Provenance columns on `avibase_id`

`avibase_id_source` records how each id was determined:

- `existing_ok` - a well-formed id that validates against the current AviList/eBird taxonomy (re-checked on every run, so a taxonomy refresh surfaces automatically)
- `ebird_single`, `ebird_group_exact`, `sci_name_ebird`, `sci_name_avilist`, `inherited_from_merge_target` - auto-derived this run (only ever seen for a *new*, not-yet-resolved row - once resolved these collapse into `existing_ok` on the next run, since the value then validates)
- `manual` - **sticky**: a human deliberately chose this id (typically a sensu-lato concept the automated cascade can't derive on its own); never silently overwritten by a re-run
- `legacy_confirmed` - **sticky**: one of 18 SEQs reviewed by hand (2026-09) whose id is absent from current AviList/eBird by design - a pre-split "sensu lato" concept where the other member does not occur in Kenya, so no further review is needed. Matched by `SEQ`, not by exact id string, so a later manual refinement of the value stays accepted. The list lives in `update_taxonomy.py` as `CONFIRMED_LEGACY_SEQ`.

The one-time pristine snapshot of `species_base_list.csv`, taken before any of
this cleanup, is kept at `species_base_list.original.csv` for audit purposes.

## Columns this script does not derive

Two groups of columns are carried through untouched, because nothing in
AviList or eBird can produce them. Both were recovered from the last commit of
`data/species_base_list.xlsx` (git history) after the rewrite to CSV dropped
them, and both are listed in `update_taxonomy.py`'s output columns so a re-run
preserves them:

- `endemic`, `afrotropical`, `palearctic`, `waterbird` - hand-curated
  ecological flags. The website filters on these.
- `mass`, `habitat`, `habitat_density`, `migration`, `trophic_level`,
  `trophic_niche`, `primary_lifestyle`, `range_size` - AVONET traits, read by
  `F_analysis.m`. `habitat_density` and `migration` are numeric codes used as
  indices into a category list, so they must stay numeric. These are **not**
  in the website export: the site never reads them.

The four AVONET coordinate columns (`Min_Latitude`, `Max_Latitude`,
`Centroid_Latitude`, `Centroid_Longitude`) were not restored - nothing
references them. They remain available in git history if ever needed.

## Historical-lump documentation

14 SEQs are lumps where at least one member was already a distinct,
field-separable species in eBird's ~2007 taxonomy - i.e. the 1970-1984 atlas
could have told them apart but didn't. These get an auto-generated `comment`
(only when `comment` was previously blank) listing the members and their
current Kenyan eBird record counts. The list (`HISTORICAL_LUMPS` in
`update_taxonomy.py`) was derived 2026-09 by diffing the current eBird
taxonomy against an archived v1.05 (~2007) pull
(`api.ebird.org/v2/ref/taxonomy/ebird?version=1.05`, not kept as a repo file);
re-derive the same way if this needs revisiting after a future split/lump.
`SEQ 786` (Collared/Pied/Semicollared Flycatcher) was reviewed the same way
but got a hand-written comment instead, since the three-way inclusion was a
deliberate editorial call, not a mechanical lump.

## What the website consumes

`E_export_website.m` writes `export/website/sp_base.json`, which is the
contract the site reads (`npm run sync:data` there copies it in, and
`tests/dataSchema.test.js` asserts it). It is a deliberately narrower view of
`species_base_list.csv`:

- the two namings travel side by side - `common_name`/`scientific_name`/`SEQ`
  for the 1989 atlas, `avilist_*` for AviList - and neither overwrites the
  other, because the site lets the reader switch between them
- `family` is exported under that name; inside the MATLAB pipeline the same
  column is aliased `checklist_family`, which `F_analysis.m` and `plot_tree.R`
  still read
- `ebird` is a flat array of species-level codes (the `|`-joined `ebird_code`
  and `avilist_members` are dropped), so every rendered link resolves. The
  previous export listed each SEQ's full eBird mapping, including slash/spuh
  codes such as `y00820`, which have no species page
- a missing value is `null`, not `""` or `"0"`

`SEQ_old`/`SEQ_new` in `map_data.json` are built with `num2cell`, so a square
always yields a flat array. The earlier `{...}` wrap produced `[[]]` for a
square with no species; the site flattens one level, so that nested empty
array landed in the SEQ set as an array object - never matching a SEQ, which
inflated the square's count. 66 squares were affected.

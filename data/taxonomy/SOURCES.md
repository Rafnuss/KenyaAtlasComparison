# Taxonomic reference sources

Reference files used to build the `SEQ` <-> AviList crosswalk that lives directly
in `data/species_base_list.csv` (see below). Re-download with
`bash data/taxonomy/fetch_sources.sh`. **Pin these versions**: the crosswalk
columns are only meaningful against the versions listed here.

## The matching rule, in full

`SEQ` is the key for everything. One row per atlas concept, and every data
source joins to it:

| Layer | Source data | Hand-edited crosswalk | Built by |
|---|---|---|---|
| 1970-84 atlas | `data/oldatlas/A Bird Atlas of Kenya_v5.xlsx` | *none - the sheet carries `SEQ` itself* | `A_import_old_atlas.m` |
| Kenya Bird Map | `data/kbm/geojson/*` | `data/kbm/sp_kbm.xlsx` (`Ref` -> `SEQ`) | `B_import_KBM.m` |
| eBird 2009-23 | `data/eBird/ebd_KE_relOct-2023/` | `data/eBird/sp_ebird.xlsx` (taxon -> `SEQ`) | `C_import_ebird.m` |
| AviList naming | `avilist_v2025b_species.csv` | *none - derived from `sp_ebird.xlsx`* | `update_taxonomy.py` |

So there are exactly **two files to hand-edit**: `sp_kbm.xlsx` for the KBM
association and `sp_ebird.xlsx` for the eBird association. Everything else in
`species_base_list.csv` - `avilist_common_name`, `avilist_scientific_name`,
`avilist_members`, `ebird_code`, `family`, `order`, `iucn`, `birdlife_url`,
`n_avilist_species`, `avibase_id` - is **derived** and is overwritten on every
`update_taxonomy.py` run. Editing those columns by hand does not stick.

`update_taxonomy.py` derives a SEQ's modern species from the eBird taxa mapped
to it in `sp_ebird.xlsx`, and then:

- **1 modern species** -> that species' name and Avibase id (`ebird_single`)
- **2+ modern species** -> eBird's own slash/spuh taxon if one covers exactly
  that set (`ebird_group_exact`), else a hand-picked sensu-lato id kept sticky
  as `manual`, else the row goes to `avibase_id_todo.csv` for a human

### Editing `sp_ebird.xlsx` safely

The file does two jobs at once, and confusing them is the one way to break it:

1. It is the **join table for the raw EBD**, matched on `scientific_name`
   (falling back to `avibase_id`). Its rows are 2023-vintage eBird taxa.
2. It **declares membership** - which modern species an atlas concept covers.

Therefore: **never delete a row whose `GroupCount` is greater than zero.** That
column is the number of Kenyan EBD records for the taxon, so deleting such a
row silently drops real observations. To re-scope a concept, *add* rows for
the current species instead. A row whose `species_code` is a `slash`/`spuh` in
the current eBird taxonomy is used for the EBD join but ignored when deriving
members, so a historical lumped taxon and its modern constituents coexist
happily on the same `SEQ`. Set `SEQ = 0` for a taxon that belongs to no atlas
concept.

SEQ 556 (Red-rumped Swallow) is the worked example: `rerswa1`
(`Cecropis daurica`, 272 Kenyan records, now a three-way slash) stays for the
join, while `rerswa8` + `rerswa12` declare the European + African membership.

### Reviewing a change

`Z_review.m` diffs the whole pipeline against a baseline commit, joined on
`SEQ`, and writes `export/review/review.html` (searchable, with per-species
diff maps) plus `review.csv`. A taxonomy edit should change names and
associations - never a map - and that page is where you confirm it.
Re-running `C_import_ebird.m` after an edit takes ~8 s thanks to
`data/eBird/ebd_reduced.mat`, so the edit/review loop is fast.

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
| `family` | AviList English family name (from the member(s)), falling back to eBird's own family when AviList has no species-rank row - `SEQ 198` (African Swamphen, which AviList ranks as a subspecies of Purple Swamphen while eBird splits it) would otherwise be blank |
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

- `ebird_single`, `ebird_group_exact`, `sci_name_ebird`, `sci_name_avilist`, `inherited_from_merge_target` - derived fresh every run directly from the row's own members (via `sp_ebird.xlsx` -> eBird taxonomy), so this is the authoritative source whenever it finds anything
- `existing_ok` - a well-formed id that doesn't derive from the row's own members (no single eBird parent, no exact slash/spuh group) but still resolves to *some* current AviList/eBird concept, so it's accepted as-is. This is deliberately the weakest check in the cascade - it only means the id isn't dangling, not that it's the id *this* row's members would produce - so it defers to member-derivation whenever that's possible (see the 2026-09 fix below)
- `manual` - **sticky**: a human deliberately chose this id (typically a sensu-lato concept the automated cascade can't derive on its own); never silently overwritten by a re-run
- `legacy_confirmed` - **sticky**: a SEQ reviewed by hand whose id is absent from current AviList/eBird by design - a pre-split "sensu lato" concept where the other member does not occur in Kenya, so no further review is needed. Matched by `SEQ`, not by exact id string, so a later manual refinement of the value stays accepted. The list lives in `update_taxonomy.py` as `CONFIRMED_LEGACY_SEQ`.

**Bug fixed 2026-09**: `CONFIRMED_LEGACY_SEQ` used to list 18 SEQs; only 1
(`991`, Rufous Sparrow) is genuinely unresolvable. 17 were misclassified
because phase 1 checked `existing_ok` (id exists *somewhere* in the global
taxonomy) *before* ever deriving what the row's own members resolve to, so a
stale id left over from years ago - still "valid" as some unrelated taxon -
was accepted and never compared against the correct one. Example: SEQ 560
(African Rock Martin) held `avibase-47DA0258`, the id for the slash taxon
"Pale/Red-throated Crag-Martin" - a different bird - instead of
`avibase-21E0ADE5`, the id its own `ebird_code` (`rocmar5`) actually resolves
to. 31 rows total had a wrong id this way (28 single-species/small-group rows
now auto-correct; 3 real 2-species lumps with no matching eBird slash/spuh -
Grey Woodpecker, Northern White-tailed Bush Lark, Abyssinian White-eye - can't
be auto-corrected and are back in `avibase_id_todo.csv` for manual input, the
same workflow as the original cleanup). Phase 1 now always derives from
members first and only falls back to `existing_ok`'s weaker check when
derivation finds nothing - and never falls back at all for a real multi-
species lump missing a group taxon, since keeping any single-member id there
would silently drop the rest of the concept.

The remaining 17th SEQ, `556` (Red-rumped Swallow), was in the set for a
different reason: `sp_ebird.xlsx` mapped it to a single old code (`rerswa1`)
that eBird's current taxonomy classifies as the 3-way slash "European/African/
Eastern Red-rumped Swallow" - broader than the atlas concept, and not
resolvable as a single member so it fell through to a blank crosswalk
entirely. Corrected by remapping it to the two real current species (European
`Cecropis rufula` + African `Cecropis melanocrissus`, excluding the Eastern/
Daurian group), which eBird's own exact "European/African Red-rumped Swallow"
slash (`y01284`, `avibase-1B08050E`) now resolves automatically - no manual
override needed any more.

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

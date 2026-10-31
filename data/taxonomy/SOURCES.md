# Taxonomy and comparison mappings

The digitization project owns the original atlas names and Avibase concepts. Comparison owns decisions about which observations can be compared. The website consumes the generated comparison data. eBird and the Birds of Kenya 2019 checklist are not naming systems in this workflow.

## Four curated tables

| Table | One row represents | Editable fields |
|---|---|---|
| `data/species_base_list.csv` | One comparison taxon, uniquely keyed by `avibase_id` | Original display names, explicit `avilist_ids` membership, ecological flags, AVONET traits and notes |
| `data/taxonomy/atlas_mapping.csv` | One historical atlas entry | Original `atlas_seq` and `atlas_avibase_id`, comparison `avibase_id`, legacy `comparison_seq`, `external_use` and decision notes |
| `data/eBird/species_mapping.csv` | One archived eBird source taxon | Source scientific name/concept ID/code, comparison `avibase_id`, current species-page `link_code`, notes |
| `data/kbm/species_mapping.csv` | One KBM `ref` | Source names, comparison `avibase_id`, notes |

A blank comparison ID means excluded or unresolved. Keep the row and its explanation. A pipe in `avilist_ids` explicitly lists the AviList members of a comparison group; neither eBird assignments nor a broad Avibase parent expands that membership automatically. `GroupCount` in the eBird table preserves the original inventory count; zero-count rows can supply modern species links without adding observations.

`SEQ` in the comparison taxon table is retained only as a legacy numeric reference for historical plots, external studies and sorting. It is not the source atlas key: Shoebill's original project code is `0a`, while its legacy comparison code is `9023`. Read `atlas_seq` and the digitization export's `seq` as text. The unnumbered Taita Thrush (`0b`) has no 1970–1984 records and no separate comparison map: the historical Northern Olive Thrush concept already includes helleri. Its empty comparison mapping is explicit.

## Sources and generated files

The sibling repository `Digitization-of-A-Bird-Atlas-of-Kenya` supplies:

- `data/curated/atlas_taxonomy.csv`: historical concepts, book-review evidence and qualifications.
- `export/comparison/atlas_records.csv`: original observations and period fields, joined to those historical IDs.
- `export/comparison/avilist_taxa.csv`: AviList v2026 fields, exported by `scripts/04_taxonomy_reference.m` from `data/raw/taxonomy/AviList-v2026-20Oct2026-extended-version.xlsx`.

Run `python3 data/taxonomy/update_taxonomy.py` to generate `data/taxonomy/comparison_taxa.csv`. This short script only reads the explicit memberships and the pinned AviList reference. It never edits the curated files, resolves IDs from names, or uses eBird to define a concept. Missing member IDs fail at the join so a new reference release requires deliberate curation.

Exact reference matches supply the AviList name. A historical group absent from AviList displays its explicit member names under “AviList members”, with `name_source` identifying the group. It is not presented as an official AviList species. A subspecies can have an AviList scientific name but no English name; the English species name comes from its explicitly listed member. Family names now use AviList's scientific family field. IUCN is the most severe available member category; no supplied assessment stays blank. Traits remain curated, not inferred from taxonomy.

The old eBird taxonomy resolver, pending-ID table, 2025b naming reference and editable Excel crosswalks have been retired. Their history remains in Git. Archived eBird observation names/IDs stay unchanged in the direct source table; `link_code` is the separately curated current page code.

## Comparison decisions

The exclusions and historical mappings remain explicit in `atlas_mapping.csv`. Barbary Falcon, Usambiro Barbet and Kilimanjaro Honeyguide retain separate historical maps, as do the two guineafowl. Modern species-level Peregrine Falcon, d'Arnaud's Barbet and Pallid Honeyguide observations can include these separate historical taxa and cannot be partitioned without finer source identifications. Other earlier folds remain under review. Atlas Broad-billed Prion and Herring Gull are reidentified as Slender-billed Prion and Lesser Black-backed Gull. These identification decisions are not applied to external-study trends (`external_use=false`).

The 2026 migration adopts the digitization project's revised IDs for the other active historical entries, including Kenyan AviList splits. Original spelling corrections and the Little/Dwarf Bittern SEQs come from the corrected atlas CSV. The historical book-review qualifications remain in the curated notes. White-eye parent IDs are not expanded into global memberships; Kulal is represented within AviList's Kafa White-eye. KBM Ref 1103 remains unassigned pending review.

Modern Melodious Lark and European Pied Flycatcher have their own comparison taxa (legacy codes 9024 and 9025). They are not members of the named historical White-tailed Lark or Collared/Semicollared Flycatcher concepts. Historical identification uncertainty remains flagged. The unresolved three-way Ficedula slash/spuh and KBM Ref 917 are retained but unassigned because they cannot be attributed to the narrower pair.

## Rebuild and review

Run `AA_wrapper_all.m` from the Comparison root. This prepares the names, imports A–C, refits the existing effort correction in D, and exports the website JSON in E. Modern maps and effort coverage are restricted to the historical atlas squares (`g.SqN > 0`), so cells outside the website map cannot inflate species totals. The raw observations and source caches are retained. Squares without historical coverage do not calibrate the effort fit; they retain an unknown correction and are masked, apart from the existing one-period-only rules. The KBM cache is indexed by its stored source `ref`; each final map layer stores its Avibase ordering, checked before combining layers. The eBird cache can be reused after mapping edits; delete it only when replacing the source download or grid.

In `KenyaBirdTrends`, run `npm run sync:data`, `npm test` and `npm run build`. Website selection, URLs and grid membership use Avibase IDs. Numeric historical URLs remain readable. The selected panel shows both naming systems; map downloads are SVGs generated from the current taxon's data, not archived PNGs.

`Z_review.m` compares maps to the 2024 baseline by legacy SEQ and displays the current direct source mappings. Changes to names and concepts can legitimately change map assignments: inspect the per-source differences before a release. This development update does not publish a website or a GBIF dataset.

## Cross-project eBird 2026 audit

Birds-of-Kenya now pins the official eBird names/taxonomy workbook dated 27 October 2026. Its `scripts/check_taxonomy_coherence.py` verifies every emitted Comparison source link against current eBird species IDs and explicit AviList membership, allowing only the documented Swamphen and Siberian Stonechat scope overrides. Current aliases supply links without adding archived observations. The source mapping retains its archived names/IDs/codes.

The book's lark accounts on printed pp. 306–307 correct two historical concepts: SEQ 531 is the former africana complex, including present Rufous-naped and Athi Larks; SEQ 532 is hypermetra/kidepoensis, including Ugandan atlas records. The previous national 2019 crosswalk incorrectly attached Athi to Red-winged Lark. Both historical group IDs and source associations now agree across projects. National checklist species remain narrower than these historical groups.

SEQ 780 retains the original Spotted/Mediterranean concept and explicitly lists both AviList 2026 members. Its archived undifferentiated records are appropriate for that broad map. Birds-of-Kenya instead uses only the source striata-group records for its current Spotted Flycatcher species; it leaves the undifferentiated records outside species evidence. Little Swift retains its historical narrower ID while linking to the broader 2026 Little/House Swift lump.

#!/usr/bin/env bash
# Re-download the taxonomic reference files used by update_taxonomy.py.
# Run from the repository root:  bash data/taxonomy/fetch_sources.sh
set -euo pipefail
cd "$(dirname "$0")"

echo "== AviList v2025b (10 Jun 2026) =="
curl -fsSL -o AviList-v2025b-10Jun2026-extended.xlsx \
  "https://www.avilist.org/wp-content/uploads/2026/06/AviList-v2025b-10Jun2026-extended.xlsx"
curl -fsSL -o AviList_v2025_metadata_11Jun.pdf \
  "https://www.avilist.org/wp-content/uploads/2025/06/AviList_v2025_metadata_11Jun.pdf"

echo "then re-extract the species-only working copy:"
echo "  (see the one-liner in SOURCES.md, or just ask to regenerate avilist_v2025b_species.csv)"

echo
echo "MANUAL STEP -- Cloudflare blocks scripted access to the Clements site,"
echo "and this is the ONLY source that carries TAXON_CONCEPT_ID (Avibase id)"
echo "for eBird's slash/spuh/hybrid taxa, so there is no automatable alternative."
echo "Download this file in a browser and save it next to this script:"
echo "  https://www.birds.cornell.edu/clementschecklist/download/"
echo "  -> 'eBird taxonomy v2025' -> CSV file"
echo "  -> save as eBird_taxonomy_v2025-4.csv"

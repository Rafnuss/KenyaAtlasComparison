#!/usr/bin/env python3
"""Build display metadata from curated comparison taxa and pinned AviList 2026.

Run from any directory. IDs, memberships and source associations are never inferred.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
with (ROOT / 'data/species_base_list.csv').open() as f:
    taxa = list(csv.DictReader(f))
with (ROOT.parent / 'Digitization-of-A-Bird-Atlas-of-Kenya/export/comparison/avilist_taxa.csv').open() as f:
    reference = {r['avibase_id']: r for r in csv.DictReader(f)}
with (ROOT / 'data/taxonomy/atlas_mapping.csv').open() as f:
    atlas_mapping = list(csv.DictReader(f))

severity = ['NE', 'DD', 'LC', 'NT', 'VU', 'EN', 'CR', 'CR (PE)', 'CR (PEW)', 'EW', 'EX']
rows = []
for taxon in taxa:
    members = [reference[i] for i in taxon['avilist_ids'].split('|') if i]
    # Only a single-member taxon takes its AviList name; any group lists its members.
    exact = reference.get(taxon['avibase_id']) if len(members) == 1 else None
    row = dict(taxon)
    row['atlas_seq'] = '|'.join(r['atlas_seq'] for r in atlas_mapping if r['avibase_id'] == taxon['avibase_id'])
    # A historical group has no exact AviList name: label its explicit members.
    row['avilist_common_name'] = exact['common_name'] if exact and exact['common_name'] else ' / '.join(m['common_name'] for m in members)
    row['avilist_scientific_name'] = exact['scientific_name'] if exact else ' / '.join(m['scientific_name'] for m in members)
    row['avilist_sort'] = min(float(m['avilist_sort']) for m in members)
    row['avilist_members'] = ' | '.join(m['scientific_name'] for m in members)
    row['n_avilist_species'] = len(members)
    row['family'] = ' | '.join(dict.fromkeys(m['family'] for m in members))
    row['order'] = ' | '.join(dict.fromkeys(m['order'] for m in members))
    statuses = [m['iucn'] for m in members if m['iucn'] in severity]
    row['iucn'] = max(statuses, key=severity.index) if statuses else ''
    row['birdlife_url'] = members[0]['birdlife_url'] if len(members) == 1 else ''
    row['name_source'] = 'AviList v2026' if exact else 'Historical concept; AviList v2026 members'
    rows.append(row)

with (ROOT / 'data/taxonomy/comparison_taxa.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
print(f'Wrote {len(rows)} comparison taxa; curated mappings unchanged.')

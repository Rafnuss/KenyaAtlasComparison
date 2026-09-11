#!/usr/bin/env python3
"""
Add an `avibase_id` column to data/eBird/sp_ebird.xlsx: the Avibase concept id
(eBird's TAXON_CONCEPT_ID) as it actually appears in this project's archived
EBD download, for joining raw EBD records to a SEQ in C_import_ebird.m.

Why not match on scientific_name (the previous approach) or species_code?
Both drift under a taxonomy update - eBird's 2025 revision renamed ~32
scientific names (genus reassignments) and recategorized/retired ~22 of this
file's species_codes (e.g. "whimbr" now names a slash, not Eurasian Whimbrel
alone). TAXON_CONCEPT_ID is the one field designed to survive that: it is
carried on every EBD record (all categories, including slash/spuh/hybrid) and
stays stable across splits/renames for a given underlying concept.

Cascade per row (first hit wins), matching update_taxonomy.py's approach:
  1. this taxon's own (category, scientific_name) as it appears somewhere in
     the archived EBD (data/eBird/ebd_KE_relOct-2023/ebd_KE_relOct-2023.txt)
     - guaranteed to be the exact id that file's own records carry
  2. this row's species_code, looked up in the current eBird taxonomy
     (data/taxonomy/eBird_taxonomy_v2025-4.csv) - for taxa never recorded in
     Kenya (2), so there's nothing to match in the archived EBD
  3. this row's scientific_name, same current-taxonomy lookup

Nothing else in sp_ebird.xlsx is touched: common_name/scientific_name/
species_code stay exactly as they were. species_base_list.csv is the
authoritative current-taxonomy crosswalk (avilist_common_name, family,
order, ...) for each SEQ; this file's job is only the raw-EBD join.

Run from the repository root:  python3 data/eBird/add_avibase_id.py
Idempotent.
"""
import csv, os, re, shutil, sys
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SP_EBIRD = os.path.join(ROOT, "data", "eBird", "sp_ebird.xlsx")
EBD_TXT = os.path.join(ROOT, "data", "eBird", "ebd_KE_relOct-2023",
                        "ebd_KE_relOct-2023.txt")
EBIRD = os.path.join(ROOT, "data", "taxonomy", "eBird_taxonomy_v2025-4.csv")
SNAPSHOT = os.path.join(ROOT, "data", "eBird", "sp_ebird.original.xlsx")

WELL_FORMED = re.compile(r"^avibase-[0-9A-F]{8}$", re.I)


def extract_ebd_concepts(path):
    """(category, scientific_name) -> avibase_id, from the raw EBD's own header."""
    m = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        i_cat = hdr.index("CATEGORY")
        i_sci = hdr.index("SCIENTIFIC NAME")
        i_tid = hdr.index("TAXON CONCEPT ID")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            key = (parts[i_cat], parts[i_sci])
            if key not in m:
                m[key] = parts[i_tid]
    return m


def main():
    for p in (SP_EBIRD, EBD_TXT, EBIRD):
        if not os.path.exists(p):
            sys.exit("missing input: " + p +
                     ("\n(extract data/eBird/ebd_KE_relOct-2023.zip first)"
                      if p == EBD_TXT else ""))

    print("scanning archived EBD for (category, scientific_name) -> avibase id "
          "(one pass over the ~900MB file)...")
    ebd_map = extract_ebd_concepts(EBD_TXT)

    with open(EBIRD, encoding="utf-8-sig") as f:
        ebird = list(csv.DictReader(f))
    eb_by_code = {r["SPECIES_CODE"]: r for r in ebird}
    eb_by_sci = {r["SCI_NAME"]: r for r in ebird}

    cats = ["species", "issf", "slash", "spuh", "hybrid", "form",
            "intergrade", "domestic"]

    wb = openpyxl.load_workbook(SP_EBIRD, data_only=True)
    ws = wb.worksheets[0]
    hdr = [c.value for c in ws[1]]
    if "avibase_id" in hdr:
        i_avid = hdr.index("avibase_id") + 1
    else:
        i_avid = len(hdr) + 1
        ws.cell(row=1, column=i_avid, value="avibase_id")
    i_sci = hdr.index("scientific_name") + 1
    i_code = hdr.index("species_code") + 1

    stats = {"archived_ebd": 0, "current_code": 0, "current_sci": 0, "none": 0}
    unresolved = []
    for row in ws.iter_rows(min_row=2):
        sci = (row[i_sci - 1].value or "").strip()
        code = (row[i_code - 1].value or "").strip()
        found = None
        for c in cats:
            if (c, sci) in ebd_map:
                found = ebd_map[(c, sci)]
                stats["archived_ebd"] += 1
                break
        if not found and code in eb_by_code:
            found = eb_by_code[code]["TAXON_CONCEPT_ID"]
            stats["current_code"] += 1
        if not found and sci in eb_by_sci:
            found = eb_by_sci[sci]["TAXON_CONCEPT_ID"]
            stats["current_sci"] += 1
        if not found:
            stats["none"] += 1
            unresolved.append((row[0].row, sci, code))
        row[i_avid - 1].value = found or ""

    if not os.path.exists(SNAPSHOT):
        shutil.copy2(SP_EBIRD, SNAPSHOT)
    wb.save(SP_EBIRD)

    n = sum(stats.values())
    print("rows: %d" % n)
    for k in ("archived_ebd", "current_code", "current_sci", "none"):
        print("  %-14s %4d" % (k, stats[k]))
    if unresolved:
        print("\nunresolved (%d):" % len(unresolved))
        for r, sci, code in unresolved:
            print("    row %d  %s / %s" % (r, sci, code))


if __name__ == "__main__":
    main()

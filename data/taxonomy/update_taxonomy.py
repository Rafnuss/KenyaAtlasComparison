#!/usr/bin/env python3
"""
Build the SEQ <-> AviList crosswalk directly inside data/species_base_list.csv.

Replaces resolve_avibase_ids.py (Step 2) and folds in Step 3 (the crosswalk)
as extra columns on the same file, rather than a separate table: every atlas
concept (SEQ) is a lump of 1+ modern species, and a lump's members almost
always share one family/order/IUCN status, so a few pipe-delimited columns on
the existing one-row-per-SEQ file carry everything without a second file.
MATLAB reads them with strsplit(x, "|").

Phases
------
  0. absorb any ids hand-filled in avibase_id_todo.csv
  1. resolve/validate avibase_id (unchanged logic from Step 2)
  2. auto-fill `comment` for historical lumps where a member was already a
     separate species in the atlas era (blank comments only - never overwrites)
  3. enrich each row with AviList-derived columns: avilist_common_name,
     avilist_scientific_name, ebird_code, family, order, iucn, n_avilist_species

Run from the repository root:  python3 data/taxonomy/update_taxonomy.py
Idempotent - re-run any time after editing avibase_id_todo.csv or AviList/eBird
reference files are refreshed.
"""
import csv, os, re, shutil, sys, zipfile
from xml.etree import ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "data", "species_base_list.csv")
AVILIST = os.path.join(ROOT, "data", "taxonomy", "avilist_v2025b_species.csv")
EBIRD = os.path.join(ROOT, "data", "taxonomy", "eBird_taxonomy_v2025-4.csv")
SP_EBIRD = os.path.join(ROOT, "data", "eBird", "sp_ebird.xlsx")
TODO = os.path.join(ROOT, "data", "taxonomy", "avibase_id_todo.csv")
SNAPSHOT = os.path.join(ROOT, "data", "taxonomy", "species_base_list.original.csv")
RECORDS = os.path.join(ROOT, "data", "taxonomy", "ebd_KE_relOct-2023_records_per_concept.tsv")

WELL_FORMED = re.compile(r"^avibase-[0-9A-F]{8}$", re.I)

# SEQs confirmed by hand (2026-09) as legitimate Avibase "sensu lato" concepts:
# each combines a pre-split pair/group where one member does not occur in
# Kenya, so an id absent from AviList/eBird is correct and deliberate here,
# not an unresolved gap. Matched by SEQ, not by exact id string, so a later
# manual refinement of the id (e.g. a closer Avibase lookup) stays accepted.
#
# Corrected 2026-09 (post-launch bug): this set used to list 18 SEQs, but only
# 1 is genuinely unresolvable - the other 17 are single/small-group, ordinary
# current species that eBird resolves cleanly via their own members; they
# only ended up here because phase 1's old validation accepted any well-formed
# id that existed *somewhere* in the global taxonomy, without checking it
# named the concept this row is actually about. Each had a stale id left over
# from years ago (still "valid" as *some* unrelated taxon, e.g. SEQ 560's
# African Rock Martin held avibase-47DA0258, the id for the slash taxon
# "Pale/Red-throated Crag-Martin" - a different bird entirely) and was then
# rubber-stamped as "legacy_confirmed" without re-deriving it from ebird_code
# first. See the corrected phase 1 below: resolve(row) now always runs first
# and wins whenever it finds a member-derived id, and the "cur exists
# somewhere" fallback only applies when resolve() found nothing to compare
# against. SEQ 556 (Red-rumped Swallow) was in this set for a different
# reason - it had no sp_ebird.xlsx member at all - but a 2026-09 fix mapped it
# to the current European + African Red-rumped Swallow (eBird's own exact
# "European/African Red-rumped Swallow" slash, y01284/avibase-1B08050E, now
# resolves it automatically), so it no longer needs this override either.
CONFIRMED_LEGACY_SEQ = {"991"}

# Historical lumps where at least one member was already a distinct species
# in eBird's ~2007 taxonomy (i.e. separable in the field during the 1970-1984
# atlas era), identified 2026-09 by cross-referencing archived eBird taxonomy
# versions (v1.05) against these SEQs' member lists and Kenyan EBD record
# counts. Re-derive by refetching an old version via
# api.ebird.org/v2/ref/taxonomy/ebird?version=1.05 if this needs revisiting.
HISTORICAL_LUMPS = {10, 100, 311, 382, 384, 386, 459, 513, 522, 528, 566, 761,
                     913, 990}


# --------------------------------------------------------------------------
# minimal dependency-free .xlsx reader (first worksheet, values as strings)
# --------------------------------------------------------------------------
def read_xlsx(path):
    NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
                shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
        name = sorted(n for n in z.namelist()
                      if re.match(r"xl/worksheets/sheet\d+\.xml$", n))[0]
        rows = []
        for row in ET.fromstring(z.read(name)).iter(NS + "row"):
            cells = {}
            for c in row.iter(NS + "c"):
                ref = c.get("r") or ""
                col = "".join(ch for ch in ref if ch.isalpha())
                n = 0
                for ch in col:
                    n = n * 26 + (ord(ch.upper()) - 64)
                v = c.find(NS + "v")
                txt = "" if v is None or v.text is None else v.text
                if c.get("t") == "s" and txt.isdigit():
                    txt = shared[int(txt)]
                elif c.get("t") == "inlineStr":
                    is_ = c.find(NS + "is")
                    txt = "".join(t.text or "" for t in is_.iter(NS + "t")) if is_ is not None else ""
                cells[n - 1] = txt.strip()
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
    if not rows:
        return []
    hdr = rows[0]
    return [dict(zip(hdr, r + [""] * (len(hdr) - len(r)))) for r in rows[1:] if any(r)]


def load_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


# most severe first; a bracketed qualifier (e.g. "CR (PE)") ranks with its
# base code, one notch more severe than the plain code
IUCN_SEVERITY = {"EX": 0, "EW": 1, "CR": 3, "EN": 4, "VU": 5, "NT": 6, "LC": 7,
                 "DD": 8, "NE": 9}


def worst_iucn(categories):
    """most severe IUCN category among a set (AviList codes, qualifiers ok)"""
    def rank(c):
        base = c.split(" (")[0]
        bump = -0.5 if "(" in c else 0  # a qualified CR ranks just above plain CR
        return IUCN_SEVERITY.get(base, 99) + bump
    cats = [c for c in categories if c]
    return min(cats, key=rank) if cats else ""


def slash_label(names):
    """eBird-style slash label: strip the longest shared trailing token run."""
    if len(names) == 1:
        return names[0]
    toks = [n.split() for n in names]
    k = 0
    while all(len(t) > k + 1 for t in toks) and len({t[-1 - k] for t in toks}) == 1:
        k += 1
    if k == 0:
        return "/".join(names)
    suffix = " ".join(toks[0][-k:])
    prefixes = [" ".join(t[:-k]) for t in toks]
    return "/".join(prefixes) + " " + suffix


def main():
    for p in (BASE, AVILIST, EBIRD, SP_EBIRD):
        if not os.path.exists(p):
            sys.exit("missing input: " + p)

    base = load_csv(BASE)

    # ---- phase 0: absorb any ids hand-filled in the TODO file --------------
    if os.path.exists(TODO):
        by_seq0 = {(r["SEQ"] or "").strip().split(".")[0]: r for r in base}
        merged = []
        for t in load_csv(TODO):
            v = (t.get("avibase_id_TO_FILL") or "").strip()
            seq = (t.get("SEQ") or "").strip().split(".")[0]
            if not v or seq not in by_seq0:
                continue
            if not WELL_FORMED.match(v):
                print("  !! SEQ %s: '%s' is not a valid avibase id - skipped" % (seq, v))
                continue
            row = by_seq0[seq]
            if not WELL_FORMED.match((row.get("avibase_id") or "").strip()):
                row["avibase_id"] = v
                merged.append(seq)
        if merged:
            print("merged %d hand-filled id(s) from %s: %s\n"
                  % (len(merged), os.path.basename(TODO), ", ".join(merged)))

    avilist = load_csv(AVILIST)
    ebird = load_csv(EBIRD)
    sp_ebird = read_xlsx(SP_EBIRD)
    records = {}
    if os.path.exists(RECORDS):
        with open(RECORDS, encoding="utf-8-sig") as f:
            next(f)
            for line in f:
                k, n = line.rstrip("\n").split("\t")
                if n.isdigit():
                    records[k] = int(n)

    # ---- lookups -------------------------------------------------------
    eb_by_code = {r["SPECIES_CODE"]: r for r in ebird}
    eb_by_sci = {r["SCI_NAME"]: r for r in ebird if r["CATEGORY"] == "species"}
    eb_by_com = {r["PRIMARY_COM_NAME"]: r for r in ebird}
    av_by_sci = {r["Scientific_name"]: r for r in avilist}
    valid_ids = ({r["TAXON_CONCEPT_ID"] for r in ebird if r["TAXON_CONCEPT_ID"]}
                 | {r["AvibaseID"] for r in avilist if r["AvibaseID"]})

    def constituents(r):
        """slash/spuh SCI_NAME -> (set of species | None for genus spuh, genus)"""
        s = r["SCI_NAME"]
        if s.endswith(" sp."):
            return None, s.split()[0]
        if "/" in s:
            g = s.split()[0]
            out = {e if " " in e else g + " " + e for e in s.split(" ", 1)[1].split("/")}
            return out, g
        return {s}, s.split()[0]

    GROUPS = [(r,) + constituents(r) for r in ebird if r["CATEGORY"] in ("slash", "spuh")]

    def find_group(want):
        genera = {w.split()[0] for w in want}
        exact = [g for g, c, _ in GROUPS if c == want]
        sup = [g for g, c, _ in GROUPS if c and c > want]
        spuh = [g for g, c, gen in GROUPS if c is None and gen in genera]
        return (exact[0] if exact else None), sup, spuh

    # ---- SEQ -> eBird members -----------------------------------------
    members = {}
    for r in sp_ebird:
        seq = (r.get("SEQ") or "").strip()
        if not seq or seq in ("0", "0.0"):
            continue
        seq = seq.split(".")[0]
        code = (r.get("species_code") or "").strip()
        hit = (eb_by_code.get(code)
               or eb_by_sci.get((r.get("scientific_name") or "").strip())
               or eb_by_com.get((r.get("common_name") or "").strip()))
        if hit:
            members.setdefault(seq, [])
            if hit not in members[seq]:
                members[seq].append(hit)

    def parents_of(seq):
        """distinct modern species making up this atlas concept, + slash/spuh + hybrid taxa mapped to it"""
        mem = members.get(seq, [])
        out, groups, hybrids = [], [], []
        for m in mem:
            if m["CATEGORY"] in ("slash", "spuh"):
                groups.append(m)
                continue
            if m["CATEGORY"] in ("hybrid", "intergrade", "domestic"):
                hybrids.append(m)
                continue
            p = eb_by_code.get(m["REPORT_AS"]) if m["REPORT_AS"] else None
            q = p or m
            if q not in out:
                out.append(q)
        return out, groups, hybrids

    def resolve(row):
        """Returns (auto_id, auto_src, needs_group).

        needs_group is True only for a genuine multi-species lump (current
        sp_ebird membership resolves to >1 modern species) where no slash/spuh
        taxon in the current eBird taxonomy covers exactly that member set.
        It signals to the caller that no fallback is safe: any id sitting in
        the row can only ever name one member, silently dropping the rest, so
        this must be escalated to a human rather than accepted as "existing".
        """
        seq = (row["SEQ"] or "").strip().split(".")[0]
        parents, groups, hybrids = parents_of(seq)
        if len(parents) == 1:
            return parents[0]["TAXON_CONCEPT_ID"], "ebird_single", False
        if len(parents) > 1:
            want = {p["SCI_NAME"] for p in parents}
            for g in groups:
                c, _ = constituents(g)
                if c == want:
                    return g["TAXON_CONCEPT_ID"], "ebird_group_exact", False
            exact, _, _ = find_group(want)
            if exact:
                return exact["TAXON_CONCEPT_ID"], "ebird_group_exact", False
            return "", "", True
        sci = (row["scientific_name"] or "").strip()
        if sci in eb_by_sci:
            return eb_by_sci[sci]["TAXON_CONCEPT_ID"], "sci_name_ebird", False
        if sci in av_by_sci and av_by_sci[sci]["AvibaseID"]:
            return av_by_sci[sci]["AvibaseID"], "sci_name_avilist", False
        return "", "", False

    # ---- phase 1: resolve avibase_id ------------------------------------
    # "manual" and "legacy_confirmed" are human judgement calls, recorded once
    # and never silently overwritten by a later run. Everything else is
    # cheap/automatic and is re-validated against the current reference data
    # every time, so a refreshed AviList/eBird taxonomy is picked up for free.
    #
    # Fixed 2026-09 (post-launch bug): resolve(row) now always runs first (for
    # any non-sticky row) and its member-derived id wins whenever it finds
    # one. The previous order checked "is cur well-formed and does it exist
    # *somewhere* in the global taxonomy" before ever deriving what the row's
    # own members resolve to - so a stale id that happened to still be valid
    # for some unrelated taxon (e.g. a different slash/spuh, or a species that
    # used to be here before a code churned) was accepted as "existing_ok" and
    # never compared against the correct one. 28 rows across single species
    # and small lumps had this exact problem; see git history for this file
    # for the audit. The weaker "cur exists somewhere" fallback now only
    # applies once resolve() has found nothing to compare against, and never
    # for a real multi-species lump missing a group taxon (needs_group) -
    # there, keeping any single-member id would silently drop the rest of the
    # concept, so it goes to avibase_id_todo.csv for a human instead.
    stats, todo = {}, []
    for row in base:
        cur = (row.get("avibase_id") or "").strip()
        recorded = (row.get("avibase_id_source") or "").strip()
        seq0 = (row["SEQ"] or "").strip().split(".")[0]

        if seq0 in CONFIRMED_LEGACY_SEQ and WELL_FORMED.match(cur):
            src = "legacy_confirmed"
        elif recorded == "manual" and WELL_FORMED.match(cur):
            src = "manual"
        else:
            auto_id, auto_src, needs_group = resolve(row)
            if auto_id:
                row["avibase_id"] = auto_id
                src = auto_src
            elif needs_group:
                row["avibase_id"] = ""
                src = "TODO"
                todo.append(row)
            elif WELL_FORMED.match(cur) and cur in valid_ids:
                src = "existing_ok"
            elif WELL_FORMED.match(cur):
                # well-formed, not derivable from members, not currently
                # valid either: the only way to reach this state is a
                # deliberate human entry (a sensu-lato pick the automated
                # cascade can't derive on its own) - so it's manual, not a
                # thing to re-flag, and becomes sticky from here on.
                src = "manual"
            else:
                row["avibase_id"] = ""
                src = "TODO"
                todo.append(row)
        row["avibase_id_source"] = src
        stats[src] = stats.get(src, 0) + 1

    by_seq = {(r["SEQ"] or "").strip().split(".")[0]: r for r in base}
    for row in list(todo):
        tgt = (row["merged_SEQ"] or "").strip().split(".")[0]
        if tgt and tgt not in ("", "0") and tgt in by_seq and by_seq[tgt]["avibase_id"]:
            row["avibase_id"] = by_seq[tgt]["avibase_id"]
            stats["TODO"] -= 1
            row["avibase_id_source"] = "inherited_from_merge_target"
            stats["inherited_from_merge_target"] = stats.get("inherited_from_merge_target", 0) + 1
            todo.remove(row)

    seen = {}
    for row in base:
        if row["avibase_id"]:
            seen.setdefault(row["avibase_id"], []).append(row)
    dups = {v: rs for v, rs in seen.items()
            if len(rs) > 1 and len({r["merged_SEQ"] for r in rs} | {r["SEQ"] for r in rs
                                                                     if r["merged_SEQ"] == by_seq.get((r["merged_SEQ"] or "").split(".")[0], {}).get("SEQ", "")}) > 1}
    # a duplicate is expected when one row's merged_SEQ points at the other
    def is_expected_dup(rs):
        seqs = {r["SEQ"] for r in rs}
        return any((r["merged_SEQ"] or "").strip().split(".")[0] in seqs for r in rs)
    dups = {v: rs for v, rs in seen.items() if len(rs) > 1 and not is_expected_dup(rs)}

    # ---- phase 2: auto-fill comments for historical lumps (blank only) ----
    n_commented = 0
    for row in base:
        seq = (row["SEQ"] or "").strip().split(".")[0]
        if not seq.isdigit() or int(seq) not in HISTORICAL_LUMPS:
            continue
        if row["comment"].strip():
            continue  # never overwrite an existing note
        parents, _, _ = parents_of(seq)
        if len(parents) < 2:
            continue
        bits = []
        for p in parents:
            n = records.get(p["TAXON_CONCEPT_ID"], 0)
            bits.append("%s (%s, %d Kenyan eBird records 2009-2023)"
                         % (p["PRIMARY_COM_NAME"], p["SCI_NAME"], n))
        row["comment"] = (
            "Historical lump of %d species, treated as one taxon (%s) in the "
            "1970-1984 atlas: %s. At least one member was already a distinct, "
            "field-separable species during the atlas era; old-atlas records "
            "for this concept may include unseparated or misidentified "
            "individuals of either." % (len(parents), row["scientific_name"], "; ".join(bits))
        )
        n_commented += 1

    # ---- phase 3: enrich with AviList-derived crosswalk columns -----------
    family_mismatch = []
    for row in base:
        seq = (row["SEQ"] or "").strip().split(".")[0]
        parents, groups, _ = parents_of(seq)
        if not parents:
            # sp_ebird.xlsx's own species_code is stale (eBird v2025 downgraded
            # it to a slash, e.g. "whimbr" -> now the Hudsonian/Eurasian slash,
            # species-level code "whimbr5"): fall back to this row's own
            # scientific_name, same as the avibase_id cascade already does.
            sci = (row["scientific_name"] or "").strip()
            hit = eb_by_sci.get(sci)
            if hit:
                parents = [hit]
            else:
                for c in ("avilist_common_name", "avilist_scientific_name",
                          "avilist_members", "ebird_code", "family", "order",
                          "iucn", "n_avilist_species", "avilist_sort",
                          "birdlife_url"):
                    row.setdefault(c, "")
                # Even with no species-level member, a slash/spuh mapped to
                # this SEQ still places it in a family - SEQ 556 resolves only
                # to the European/African/Eastern Red-rumped Swallow slash,
                # which is unambiguously Hirundinidae. Worth filling: leaving
                # it blank makes every consumer special-case a missing family.
                fam = {g["FAMILY"].split("(", 1)[1].rstrip(")")
                       for g in groups if "(" in g["FAMILY"]}
                ordr = {g["ORDER"] for g in groups if g["ORDER"]}
                if len(fam) == 1:
                    row["family"] = fam.pop()
                if len(ordr) == 1:
                    row["order"] = ordr.pop()
                continue

        av = [av_by_sci.get(p["SCI_NAME"]) for p in parents]
        av = [a for a in av if a]

        # display name: prefer the eBird slash's own name if one covers this
        # exact member set (matches avibase_id resolution), else generate one
        want = {p["SCI_NAME"] for p in parents}
        exact, _, _ = find_group(want) if len(parents) > 1 else (None, None, None)
        if exact:
            row["avilist_common_name"] = exact["PRIMARY_COM_NAME"]
        elif av:
            row["avilist_common_name"] = slash_label([a["English_name_AviList"] for a in av])
        else:
            row["avilist_common_name"] = slash_label([p["PRIMARY_COM_NAME"] for p in parents])

        # avilist_members: every modern species making up this concept,
        # "|"-joined. avilist_scientific_name is the single displayable form
        # that mirrors it: the binomial itself for one species, else the
        # shared genus + "sp." (e.g. "Ficedula sp.") when the members agree
        # on a genus, else the binomials joined with "/". Keeping the pair
        # here rather than in A_import_old_atlas.m means the website and the
        # MATLAB pipeline get the same string without computing it twice.
        member_sci = [p["SCI_NAME"] for p in parents]
        row["avilist_members"] = "|".join(member_sci)
        if len(member_sci) == 1:
            row["avilist_scientific_name"] = member_sci[0]
        else:
            genera = {m.split()[0] for m in member_sci}
            row["avilist_scientific_name"] = (genera.pop() + " sp." if len(genera) == 1
                                              else "/".join(member_sci))
        row["ebird_code"] = "|".join(p["SPECIES_CODE"] for p in parents)
        row["n_avilist_species"] = str(len(parents))

        families = {a["Family_English_name"] for a in av if a["Family_English_name"]}
        orders = {a["Order"] for a in av if a["Order"]}
        # Fall back to eBird's own family/order when AviList has no
        # species-rank row for this taxon. Two concepts need it: SEQ 198
        # (African Swamphen, which AviList ranks as a subspecies of Purple
        # Swamphen while eBird splits it - one of the documented ~33
        # disagreements) and SEQ 556. Leaving them blank made every consumer
        # handle a missing family: F_analysis.m groups by it and builds a
        # filename from the group, which a blank silently breaks.
        if not families:
            # eBird spells it "Rallidae (Rails, Gallinules, and Coots)";
            # AviList's Family_English_name is just the parenthetical part.
            families = {p["FAMILY"].split("(", 1)[1].rstrip(")")
                        for p in parents if "(" in p["FAMILY"]}
        if not orders:
            orders = {p["ORDER"] for p in parents if p["ORDER"]}
        row["family"] = sorted(families)[0] if len(families) == 1 else "|".join(sorted(families))
        row["order"] = sorted(orders)[0] if len(orders) == 1 else "|".join(sorted(orders))
        if len(families) > 1 or len(orders) > 1:
            family_mismatch.append((seq, row["common_name"], families, orders))

        # IUCN: for a lump, report the most severe status among members (the
        # usual conservation-reporting convention - e.g. SEQ 753 Bar-throated
        # Apalis includes CR Apalis fuscigularis, so the concept is CR, not
        # blank, since that member alone would otherwise silently vanish from
        # any downstream by-threat-category analysis).
        cats = [a["IUCN_Red_List_Category"] for a in av if a["IUCN_Red_List_Category"]]
        row["iucn"] = worst_iucn(cats)

        # avilist_sort: AviList's own linear sequence, for sorting the website
        # species list "by taxonomy" under the AviList naming. A lump sorts at
        # its earliest member's position.
        seqs = [int(a["Sequence"]) for a in av if a["Sequence"]]
        row["avilist_sort"] = str(min(seqs)) if seqs else ""

        # birdlife_url: only for an unambiguous single species, same rule as
        # iucn - a lump has no single defensible factsheet to link to.
        row["birdlife_url"] = av[0]["BirdLife_DataZone_URL"] if len(av) == 1 else ""

    # ---- write ---------------------------------------------------------
    if not os.path.exists(SNAPSHOT):
        shutil.copy2(BASE, SNAPSHOT)
    cols = ["SEQ", "common_name", "scientific_name", "merged_SEQ", "ADU",
            "avibase_id", "avibase_id_source",
            "avilist_common_name", "avilist_scientific_name", "avilist_sort",
            "avilist_members", "ebird_code", "family", "order", "iucn",
            "birdlife_url", "n_avilist_species",
            # Not touched by this script, and not derivable from AviList.
            # The ecological flags are hand-curated; the traits come from
            # AVONET. Both were recovered from the last commit of
            # data/species_base_list.xlsx (git history) after the CSV rewrite
            # dropped them - F_analysis.m reads the traits. They stay out of
            # the website export, which only needs the taxonomy.
            "endemic", "afrotropical", "palearctic", "waterbird",
            "mass", "habitat", "habitat_density", "migration",
            "trophic_level", "trophic_niche", "primary_lifestyle",
            "range_size",
            "flag", "comment"]
    with open(BASE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in base:
            w.writerow(row)

    with open(TODO, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["SEQ", "common_name", "scientific_name", "merged_SEQ",
                    "n_members", "member_species", "member_avibase_ids",
                    "candidate_superset", "candidate_genus_spuh",
                    "avibase_search", "why", "avibase_id_TO_FILL"])
        for row in todo:
            seq = (row["SEQ"] or "").strip().split(".")[0]
            parents, groups, hybrids = parents_of(seq)
            names = [p["SCI_NAME"] for p in parents]
            ids = [p["TAXON_CONCEPT_ID"] for p in parents]
            cand_sup = cand_spuh = ""
            if names:
                _, sup, spuh = find_group(set(names))
                cand_sup = " | ".join("%s = %s (%s)" % (g["TAXON_CONCEPT_ID"], g["PRIMARY_COM_NAME"], g["SCI_NAME"]) for g in sup)
                cand_spuh = " | ".join("%s = %s" % (g["TAXON_CONCEPT_ID"], g["PRIMARY_COM_NAME"]) for g in spuh)
            why = ("no eBird taxon mapped to this SEQ" if not parents else
                   "%d modern species, no exact slash/spuh in eBird" % len(parents))
            if hybrids:
                why += "; also mapped to hybrid " + ",".join(h["SPECIES_CODE"] for h in hybrids)
            w.writerow([row["SEQ"], row["common_name"], row["scientific_name"],
                        row["merged_SEQ"], len(parents), " | ".join(names), " | ".join(ids),
                        cand_sup, cand_spuh,
                        "https://avibase.bsc-eoc.org/search.jsp?qstr=" + (row["scientific_name"] or "").replace(" ", "+"),
                        why, ""])

    # ---- report --------------------------------------------------------
    print("rows: %d" % len(base))
    for k in sorted(stats, key=lambda k: -stats[k]):
        print("  %-24s %4d" % (k, stats[k]))
    print("\nstill TODO           : %d  -> %s" % (len(todo), os.path.relpath(TODO, ROOT)))
    print("legacy-confirmed ids  : %d (accepted sensu lato concepts, not flagged)"
          % stats.get("legacy_confirmed", 0))
    if dups:
        print("\nunexpected duplicate avibase_id (%d):" % len(dups))
        for v, rs in dups.items():
            print("    %s  <- %s" % (v, ", ".join("SEQ %s %s" % (r["SEQ"], r["common_name"]) for r in rs)))
    print("\nauto-commented historical lumps this run: %d" % n_commented)
    if family_mismatch:
        print("\nfamily/order disagreement among members (%d) - check manually:" % len(family_mismatch))
        for seq, nm, fam, ordr in family_mismatch:
            print("    SEQ %s %s: families=%s orders=%s" % (seq, nm, fam, ordr))


if __name__ == "__main__":
    main()

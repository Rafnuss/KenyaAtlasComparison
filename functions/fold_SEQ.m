function seq = fold_SEQ(seq)
%FOLD_SEQ  Re-point an external dataset's atlas SEQs onto current concepts.
%
%   Several files under data/ were keyed to the 1989 atlas by hand (SABAP1-2,
%   burns2021, ...) and still carry SEQs that species_base_list.csv has since
%   retired through merged_SEQ. A_import_old_atlas.m folds the atlas records
%   themselves onto the surviving concept, so a retired SEQ matches nothing
%   downstream and an unfolded join drops that row in silence. Call this on
%   the external table's SEQ column before joining; the usual
%   `t = t(t.SEQ>0, :)` filter then discards whatever should not join.
%
%   Only *taxonomic* folds are followed. The other folds are identification
%   decisions (see each row's comment in species_base_list.csv): there the
%   external row refers to a genuinely different, still-valid species that
%   simply has no Kenyan counterpart, so following the fold would attribute
%   two foreign trends to one Kenyan concept - burns2021, for instance,
%   reports European Herring Gull (SEQ 299) and Lesser Black-backed Gull
%   (SEQ 298) separately, while in Kenya the Herring Gull records were
%   re-identified as the latter. Those SEQs become 0, i.e. dropped, as do the
%   merged_SEQ == 0 rows (historical records rejected as misidentifications).

% Classification of the merged_SEQ > 0 rows, read off their comment column.
TAXONOMIC = [144 498 502 506 691]';      % comment starts "Taxonomic: ..."
IDENTIFICATION = [299 348]';             % comment starts "NOT taxonomic ..."

sp_base = readtable("data/species_base_list.csv", TextType="string");

% Guard against a fold being added to the CSV but not classified here.
unclassified = setdiff(sp_base.SEQ(sp_base.merged_SEQ > 0), [TAXONOMIC; IDENTIFICATION]);
assert(isempty(unclassified), "fold_SEQ: SEQ %s carries a merged_SEQ but is " + ...
    "classified neither taxonomic nor identification - read its comment in " + ...
    "data/species_base_list.csv and add it to the matching list in this file.", ...
    mat2str(unclassified'));

[found, loc] = ismember(TAXONOMIC, sp_base.SEQ);
assert(all(found), "fold_SEQ: a SEQ listed as TAXONOMIC is absent from species_base_list.csv")
target = sp_base.merged_SEQ(loc);

seq(ismember(seq, [sp_base.SEQ(sp_base.merged_SEQ == 0); IDENTIFICATION])) = 0;

[tf, loc] = ismember(seq, TAXONOMIC);
seq(tf) = target(loc(tf));
end

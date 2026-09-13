%% Z_review.m - review the taxonomy update against a known-good baseline
%
% One question, asked once per SEQ: did this update change anything other
% than the species' *name* and its *eBird / KBM association*?
%
% SEQ is the key for everything - the 1970-1984 atlas sheet, the KBM
% crosswalk (data/kbm/sp_kbm.xlsx), the eBird crosswalk
% (data/eBird/sp_ebird.xlsx) and the species list
% (data/species_base_list.csv) all join on it - so this review joins on it
% too, never on row position. That matters: the baseline's sp_base is not
% even sorted by SEQ, and both lists have gained and lost rows since.
%
% Outputs, both rewritten from scratch on every run:
%   export/review/review.html  browsable - search, filter to changes only,
%                              and a per-SEQ diff map for each of the three
%                              layers (old atlas / eBird / KBM)
%   export/review/review.csv   the same table, for sorting in Excel
%
% Run after A_import_old_atlas, B_import_KBM and C_import_ebird.

% 2024-08-30 "Add data": the last commit where all three map layers were
% generated together from one species list, so the only fully self-consistent
% state to compare against. (The 2026-09 AviList commit rebuilt oldatlas.mat
% and ebirdatlas.mat but not kbmatlas.mat, which is what left the KBM layer
% keyed to a stale species list.)
BASELINE = "73f27c1";

OUT = "export/review";
if ~isfolder(OUT), mkdir(OUT); end

%% Fetch the baseline maps out of git (once; delete data/baseline to refresh)
bdir = "data/baseline";
if ~isfolder(bdir), mkdir(bdir); end
for f = ["oldatlas", "ebirdatlas", "kbmatlas"]
    p = fullfile(bdir, f + ".mat");
    if ~isfile(p)
        [st, out] = system(sprintf('git show %s:data/%s.mat > "%s"', BASELINE, f, p));
        assert(st == 0, "cannot extract %s from %s: %s", f, BASELINE, out);
    end
end

%% Load both sides
B  = load(fullfile(bdir, "oldatlas.mat"));      % map_old, sp_base (2024)
Be = load(fullfile(bdir, "ebirdatlas.mat"));    % map_ebird
Bk = load(fullfile(bdir, "kbmatlas.mat"));      % map_kbm
N  = load("data/oldatlas.mat");                 % map_old, sp_base (current)
Ne = load("data/ebirdatlas.mat");
Nk = load("data/kbmatlas.mat");
load("data/grid", "g")

% Every layer on a side must be indexed by that side's own species list.
assert(size(B.map_old, 3) == height(B.sp_base) && size(Be.map_ebird, 3) == height(B.sp_base) ...
    && size(Bk.map_kbm, 3) == height(B.sp_base), "baseline layers disagree on species count")
assert(size(N.map_old, 3) == height(N.sp_base) && size(Ne.map_ebird, 3) == height(N.sp_base) ...
    && size(Nk.map_kbm, 3) == height(N.sp_base), ...
    "current layers disagree on species count - re-run A, B and C together")
assert(isequal(size(B.map_old, [1 2]), size(N.map_old, [1 2])), "grid changed size")

%% Side tables used only for display
sp   = readtable("data/species_base_list.csv", 'TextType', 'string');
kbm  = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');
ebx  = readtable("data/eBird/sp_ebird.xlsx", 'TextType', 'string');

%% Compare, joined on SEQ
allseq = union(B.sp_base.SEQ, N.sp_base.SEQ);
allseq = allseq(~isnan(allseq));
n = numel(allseq);

[~, ib] = ismember(allseq, B.sp_base.SEQ);   % 0 where the SEQ is new
[~, in] = ismember(allseq, N.sp_base.SEQ);   % 0 where the SEQ was dropped

layers = ["old", "ebird", "kbm"];
Bmaps  = {B.map_old, Be.map_ebird, Bk.map_kbm};
Nmaps  = {N.map_old, Ne.map_ebird, Nk.map_kbm};

cnt = struct();
for L = 1:3
    cnt.(layers(L)) = zeros(n, 4);   % [base now added removed]
end
hex = strings(n, 6);                 % [old_b old_n ebird_b ebird_n kbm_b kbm_n]

for i = 1:n
    for L = 1:3
        a = getmask(Bmaps{L}, ib(i));
        b = getmask(Nmaps{L}, in(i));
        cnt.(layers(L))(i, :) = [sum(a) sum(b) sum(b & ~a) sum(a & ~b)];
        hex(i, 2 * L - 1) = tohex(a);
        hex(i, 2 * L)     = tohex(b);
    end
end

%% Assemble the review table
T = table();
T.SEQ = allseq;

T.atlas_name     = pick(N.sp_base.common_name, in, "");
T.atlas_sci      = pick(N.sp_base.scientific_name, in, "");
T.base_name      = pick(B.sp_base.common_name, ib, "");
T.base_sci       = pick(B.sp_base.scientific_name, ib, "");
% The two taxonomies the baseline carried, so a rename is reviewable against
% what the species used to be called rather than only against the atlas name.
T.base_clements  = pick(B.sp_base.clements_common_name, ib, "");
T.base_checklist = pick(B.sp_base.checklist_common_name, ib, "");

[~, is] = ismember(allseq, sp.SEQ);
T.avilist_name   = pick(sp.avilist_common_name, is, "");
T.avilist_sci    = pick(sp.avilist_scientific_name, is, "");
T.n_members      = pick(sp.n_avilist_species, is, "");
T.ebird_code     = pick(sp.ebird_code, is, "");
T.avibase_id     = pick(sp.avibase_id, is, "");
T.id_source      = pick(sp.avibase_id_source, is, "");
T.merged_SEQ     = pick(sp.merged_SEQ, is, "");
T.comment        = pick(sp.comment, is, "");

% eBird taxa and KBM references currently pointed at each SEQ. Both are the
% hand-edited crosswalks, so showing them here is what makes an association
% mistake reviewable without opening the spreadsheets.
T.ebird_mapped = joinby(allseq, ebx.SEQ, ebx.species_code + " (" + ebx.common_name + ")");
T.kbm_mapped   = joinby(allseq, kbm.SEQ, string(kbm.Ref) + " (" + kbm.Common_species + " " + kbm.Common_group + ")");

for L = 1:3
    c = cnt.(layers(L));
    T.(layers(L) + "_base")    = c(:, 1);
    T.(layers(L) + "_now")     = c(:, 2);
    T.(layers(L) + "_added")   = c(:, 3);
    T.(layers(L) + "_removed") = c(:, 4);
end

changed = cnt.old(:, 3) + cnt.old(:, 4) + cnt.ebird(:, 3) + cnt.ebird(:, 4) ...
        + cnt.kbm(:, 3) + cnt.kbm(:, 4) > 0;

status = repmat("ok", n, 1);
status(changed) = "map_changed";
status(ib == 0) = "new_seq";
status(in == 0) = "dropped_seq";
T.status = status;

% The atlas name is historical: it describes what the 1970-1984 observers
% recorded, so no taxonomy update has any business changing it. Flagged
% separately from the map diff because it means a different kind of mistake.
T.atlas_renamed = ib > 0 & in > 0 & T.atlas_name ~= T.base_name;

T = sortrows(T, "SEQ");

%% Write the CSV
writetable(T, fullfile(OUT, "review.csv"), 'QuoteStrings', true);

%% Write the HTML
% One JSON object per SEQ, carrying its own map bitmasks, so the page needs
% no column-index bookkeeping. (jsonencode flattens a 2-D cell array, so
% table2cell is not an option here.)
R = T;
for v = string(R.Properties.VariableNames)
    if isstring(R.(v))
        c = R.(v); c(ismissing(c)) = ""; R.(v) = c;
    end
end
R.old_b = hex(:, 1); R.old_n = hex(:, 2);
R.eb_b  = hex(:, 3); R.eb_n  = hex(:, 4);
R.kbm_b = hex(:, 5); R.kbm_n = hex(:, 6);

data = struct();
data.grid      = struct('nlat', numel(g.lat), 'nlon', numel(g.lon), ...
                        'valid', tohex(reshape(g.SqN > 0, [], 1)));
data.baseline  = BASELINE;
data.generated = string(datetime("now", 'Format', 'yyyy-MM-dd HH:mm'));
data.rows      = table2struct(R);

tpl  = fileread("review/template.html");
html = string(strrep(tpl, "/*__DATA__*/", jsonencode(data)));
fid = fopen(fullfile(OUT, "review.html"), 'w');
fprintf(fid, '%s', html);
fclose(fid);

% Same page, minus the document wrapper, for publishing as a shareable
% Artifact: that publisher supplies its own <head> and expects the file to
% start at the page content, so a file carrying its own <!doctype>/<html>/
% <head>/<body> is rejected. <title> and <style> stay - they are what names
% the page and themes it.
art = erase(html, ["<!doctype html>", "<html lang=""en"">", "<head>", ...
                   "</head>", "<body>", "</body>", "</html>"]);
art = regexprep(art, '<meta[^>]*>', '');
fid = fopen(fullfile(OUT, "review_artifact.html"), 'w');
fprintf(fid, '%s', art);
fclose(fid);

%% Report
fprintf('\n=== review vs %s ===\n', BASELINE);
fprintf('SEQ compared        : %d  (baseline %d, current %d)\n', n, height(B.sp_base), height(N.sp_base));
fprintf('  new SEQ           : %d\n', sum(ib == 0));
fprintf('  dropped SEQ       : %d\n', sum(in == 0));
fprintf('  atlas name changed: %d   <- should be 0\n', sum(T.atlas_renamed));
for L = 1:3
    c = cnt.(layers(L));
    both = ib > 0 & in > 0;
    ch = both & (c(:, 3) + c(:, 4) > 0);
    fprintf('  %-6s map changed : %4d SEQ  (+%d / -%d squares)\n', ...
        layers(L), sum(ch), sum(c(ch, 3)), sum(c(ch, 4)));
end
fprintf('\nwrote %s, review.csv and review_artifact.html\n', fullfile(OUT, "review.html"));

%% ---------------------------------------------------------------- helpers
function m = getmask(map, i)
% Species i of a layer as a flat logical column; all-false when absent.
if i > 0
    m = reshape(map(:, :, i), [], 1) > 0;
else
    m = false(size(map, 1) * size(map, 2), 1);
end
end

function h = tohex(v)
% Logical vector -> hex string, 4 cells per character, for the HTML page.
v = double(v(:))';
v(end + 1:ceil(numel(v) / 4) * 4) = 0;
nib = reshape(v, 4, [])' * [8; 4; 2; 1];
h = string(reshape(lower(dec2hex(nib, 1))', 1, []));
end

function out = pick(col, idx, dflt)
% col(idx), with dflt wherever idx is 0. Keeps strings as strings.
out = repmat(string(dflt), numel(idx), 1);
out(idx > 0) = string(col(idx(idx > 0)));
out(ismissing(out)) = string(dflt);
end

function out = joinby(keys, srckeys, srctext)
% All srctext whose srckey equals each key, " | "-joined.
out = strings(numel(keys), 1);
for i = 1:numel(keys)
    hit = srckeys == keys(i);
    out(i) = strjoin(srctext(hit), " | ");
end
end

% Load necessary grid and atlas data
load('data/grid')
load('data/oldatlas')

%% Read and filter eBird data
%
% Reading, filtering and gridding the 864 MB EBD text file takes minutes and
% depends only on the raw download and the grid - never on the taxonomy. It
% is therefore cached, so that re-running this script after editing
% data/eBird/sp_ebird.xlsx (the only input a taxonomy review actually
% changes) takes seconds instead of minutes.
%
% Delete data/eBird/ebd_reduced.mat to force a re-read - needed only after a
% new EBD release or a change to the grid.
CACHE = "data/eBird/ebd_reduced.mat";

if isfile(CACHE)
    load(CACHE, "ebd", "coverage_ebird")
else
    % Read the eBird data from a text file
    ebd0 = readtable("data/eBird/ebd_KE_relOct-2023/ebd_KE_relOct-2023.txt", 'TextType', 'string');
    ebd = ebd0;

    % Filter data: Keep observations from 2009 to 2023
    ebd = ebd(year(ebd.OBSERVATIONDATE) >= 2009 & year(ebd.OBSERVATIONDATE) <= 2023, :);

    % Filter out checklists that cover more than 30 km
    numel(unique(ebd.SAMPLINGEVENTIDENTIFIER(ebd.EFFORTDISTANCEKM > 30)))
    ebd(ebd.EFFORTDISTANCEKM > 30, :) = [];

    % Map eBird observations to the nearest grid cells
    [~, id_lat] = min((g.lat - ebd.LATITUDE) .^ 2, [], 2);
    [~, id_lon] = min((g.lon - ebd.LONGITUDE) .^ 2, [], 2);
    ebd.idg = sub2ind(size(g.LAT), id_lat, id_lon);

    % Compute coverage map: summarize by checklist, then by grid cell
    ebd_checklist = groupsummary(ebd, {'SAMPLINGEVENTIDENTIFIER', 'ALLSPECIESREPORTED', 'DURATIONMINUTES', 'EFFORTDISTANCEKM', 'PROTOCOLTYPE', 'NUMBEROBSERVERS', 'idg'});
    ebd_grid = groupsummary(ebd_checklist, "idg", "sum", "DURATIONMINUTES");

    coverage_ebird = nan(size(g.LAT));
    coverage_ebird(ebd_grid.idg) = ebd_grid.sum_DURATIONMINUTES / 60; % Convert minutes to hours
    coverage_ebird(isnan(coverage_ebird)) = 0;

    % Select and simplify eBird data for further processing
    % avibase_id (eBird's TAXON CONCEPT ID) is carried through for taxonomy
    % matching below: unlike scientific_name/common_name it is stable across
    % eBird taxonomy versions, so the join does not break when eBird revises
    % names or splits/lumps species (see data/eBird/add_avibase_id.py).
    ebd = table(ebd.LATITUDE, ebd.LONGITUDE, ebd.COMMONNAME, ebd.SCIENTIFICNAME, ebd.CATEGORY, ebd.TAXONCONCEPTID, ebd.idg, 'VariableNames', {'lat', 'lon', 'common_name', 'scientific_name', 'category', 'avibase_id', 'idg'});
    ebd = unique(ebd, "sorted");

    save(CACHE, "ebd", "coverage_ebird")
end

%% Taxonomy Matching

% Load the species taxonomy data for eBird
sp_ebird = readtable('data/eBird/sp_ebird.xlsx', 'TextType', 'string');

% Match species in eBird data with those in the existing species base.
% Primary key: scientific_name text, exactly as before - this is the RIGHT
% granularity for "issf" (identifiable-subspecies-group) records, since
% several distinct issf concepts (each with its own avibase_id) commonly
% share one species-level scientific_name, and sp_ebird.xlsx only tracks one
% id per species. avibase_id alone is too fine-grained here: matching on it
% first left 12,564 species/issf records unmatched (2026-09 testing).
% Fallback: avibase_id, for any record the name match misses - this is what
% makes the join robust to a *future* re-run on a newer EBD release, where
% eBird may have renamed a genus (this sp_ebird.xlsx keeps its 2023-vintage
% names) but the underlying concept id is unchanged.
[Lia_name, Locb_name] = ismember(ebd.scientific_name, sp_ebird.scientific_name);
[Lia_id, Locb_id] = ismember(ebd.avibase_id, sp_ebird.avibase_id);
Lia = Lia_name | Lia_id;
Locb = Locb_name;
Locb(~Lia_name & Lia_id) = Locb_id(~Lia_name & Lia_id);

% Ensure that all species are matched
tmp2 = unique(ebd((ebd.category == "species" | ebd.category == "issf") & ~Lia, ["common_name", "scientific_name"]));
assert(height(tmp2) == 0)

% Report species that were not matched and are categorized as "slash"
disp("Species not kept (eBird): We should only have slash ")
unique(ebd(ebd.category == "slash" & ~Locb, ["common_name", "scientific_name"]))

% Keep only matched species
ebd = ebd(Lia, :);

% Add SEQ number to the eBird data based on the species base
ebd.SEQ = sp_ebird.SEQ(Locb(Lia));

% Check for species with SEQ number 0
unique(ebd.common_name(ebd.SEQ == 0))

%% Spatial Grid Processing

% Map species observations to grid cells
[~, id_sp] = ismember(ebd.SEQ, sp_base.SEQ);
[~, id_lat] = min((g.lat - ebd.lat) .^ 2, [], 2);
[~, id_lon] = min((g.lon - ebd.lon) .^ 2, [], 2);

% Initialize an empty map for eBird data
map_ebird = false(size(map_old));

% Fill the map with species presence data
id = sub2ind(size(map_ebird), id_lat(id_sp > 0), id_lon(id_sp > 0), id_sp(id_sp > 0));
map_ebird(id) = true;

%% Save the eBird Atlas Data

% Save the processed eBird map and coverage data to a MAT-file. seq_map is
% the species key for the third dimension - see the note in B_import_KBM.m.
seq_map = sp_base.SEQ;
save('data/ebirdatlas.mat', "map_ebird", "coverage_ebird", "seq_map")
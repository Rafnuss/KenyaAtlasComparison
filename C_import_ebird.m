% Load necessary grid and atlas data
load('data/grid')
load('data/oldatlas')

%% Read and filter eBird data
%
% Reading, filtering and gridding the 864 MB EBD text file takes minutes and
% depends only on the raw download and the grid - never on the taxonomy. It
% is therefore cached, so that re-running this script after editing
% data/eBird/species_mapping.csv (the only input a taxonomy review actually
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

%% Join the raw eBird names to the curated comparison IDs
sp_ebird = readtable('data/eBird/species_mapping.csv', TextType='string');
% The download uses species-level names for identifiable subspecies groups.
% Match those source names directly; use source concept IDs for remaining taxa.
[matched_name, idx] = ismember(ebd.scientific_name, sp_ebird.scientific_name);
[matched_id, idx_id] = ismember(ebd.avibase_id, sp_ebird.source_avibase_id);
idx(~matched_name & matched_id) = idx_id(~matched_name & matched_id);
matched = matched_name | matched_id;
assert(~any((ebd.category == "species" | ebd.category == "issf") & ~matched));
ebd = ebd(matched, :);
ebd.comparison_id = sp_ebird.avibase_id(idx(matched));

%% Spatial Grid Processing

% Map species observations to grid cells
[~, id_sp] = ismember(ebd.comparison_id, sp_base.avibase_id);
[~, id_lat] = min((g.lat - ebd.lat) .^ 2, [], 2);
[~, id_lon] = min((g.lon - ebd.lon) .^ 2, [], 2);

% Initialize an empty map for eBird data
map_ebird = false(size(map_old));

% Fill the map with species presence data
id = sub2ind(size(map_ebird), id_lat(id_sp > 0), id_lon(id_sp > 0), id_sp(id_sp > 0));
map_ebird(id) = true;

% Exclude rectangular-grid cells outside the exported atlas squares.
map_ebird = map_ebird & (g.SqN > 0);
coverage_ebird(g.SqN == 0) = 0;

%% Save the eBird Atlas Data

% Save the processed eBird map and coverage data to a MAT-file. seq_map is
% the species key for the third dimension - see the note in B_import_KBM.m.
seq_map = sp_base.SEQ;
avibase_map = sp_base.avibase_id;
save('data/ebirdatlas.mat', "map_ebird", "coverage_ebird", "seq_map", "avibase_map")
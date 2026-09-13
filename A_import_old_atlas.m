%% Load and Process Grid Data

% Read the grid data generated from the R code (stored in a geojson file)
grid = jsondecode(fileread('data/oldatlas/grid.geojson'));

% Retrieve the coordinates of the center of each cell in the grid
coord = nan(2, numel(grid.features));  % Initialize a matrix to store coordinates
for i_f = 1:numel(grid.features)
    % Calculate the mean of the coordinates for each grid cell and round to one decimal
    coord(:, i_f) = round(squeeze(mean(grid.features(i_f).geometry.coordinates)), 1);
end

% Define the grid resolution and longitude/latitude values
g.res = min(diff(sort(unique(coord(1,:)))));  % Resolution of the grid
g.lon = round(min(coord(1,:)):g.res:max(coord(1,:)), 1);  % Longitude values
g.lat = round(min(coord(2,:)):g.res:max(coord(2,:)), 1);  % Latitude values
[g.LON, g.LAT] = meshgrid(g.lon, g.lat);  % Create a meshgrid for longitude and latitude

% Initialize matrices to store square labels (SqL), square numbers (SqN), and coverage
g.SqL = strings(numel(g.lat), numel(g.lon));  % Square labels
g.SqN = zeros(numel(g.lat), numel(g.lon));    % Square numbers
coverage_old = strings(numel(g.lat), numel(g.lon));  % Coverage data

% Populate the matrices with data from the grid
for i_f = 1:numel(grid.features)
    % Find the indices for latitude and longitude
    id_lat = coord(2, i_f) == g.lat;
    id_lon = g.lon == coord(1, i_f);
    
    % Ensure that only one matching index is found for both latitude and longitude
    assert(sum(id_lat) == 1 & sum(id_lon) == 1)
    
    % Assign the corresponding SqL, SqN, and coverage values
    g.SqL(id_lat, id_lon) = grid.features(i_f).properties.SqL;
    g.SqN(id_lat, id_lon) = grid.features(i_f).properties.SqN;
    coverage_old(id_lat, id_lon) = grid.features(i_f).properties.coverage;
end

% Convert coverage data from strings to numerical values
cov = nan(size(coverage_old));
cov(coverage_old == "0") = 0;
cov(coverage_old == "1-10") = 1;
cov(coverage_old == "11-30") = 2;
cov(coverage_old == "31-50") = 3;
cov(coverage_old == "51-75") = 4;
cov(coverage_old == "75-100") = 5;

% Display the coverage map
figure; 
imagesc(cov); 
set(gca, 'ydir', 'normal')  % Ensure the y-axis is in the correct direction

%% Import and Process Atlas Data

% Read the old atlas data from an Excel file
old_atlas = readtable("data/oldatlas/A Bird Atlas of Kenya_v5.xlsx", 'TextType', 'string');

% Retain only the columns of interest: SEQ, SqN, SqL, pre_1970, and x1970_1984
old_atlas = old_atlas(:, ["SEQ", "SqN", "SqL", "pre_1970", "x1970_1984"]);

% Remove records with no data prior to 1984 (pre_1970 and x1970_1984 columns)
old_atlas(ismissing(old_atlas.pre_1970) & ismissing(old_atlas.x1970_1984), :) = [];

% Retain only records with data from the main period (1970-1984)
old_atlas(ismissing(old_atlas.x1970_1984), :) = [];

% Match grid squares from the old atlas to the processed grid
[~, old_atlas.idg] = ismember(string(old_atlas.SqN) + old_atlas.SqL, string(g.SqN(:)) + g.SqL(:));

%% Create Dataset for Comparison with New Atlas

% Load the species base list. Since 2026-09 this carries the AviList/eBird
% taxonomy crosswalk built by data/taxonomy/update_taxonomy.py: avibase_id,
% avilist_common_name, avilist_scientific_name, avilist_sort, avilist_members,
% ebird_code, family, order, iucn, birdlife_url, n_avilist_species - plus the
% hand-curated ecological flags and the AVONET traits, which that script
% carries through untouched (see data/taxonomy/SOURCES.md).
sp_base = readtable("data/species_base_list.csv", 'TextType', 'string');

% Remove records with merged_SEQ == 0: these are rejected/misidentified
% historical records (see the comment column for the specific rationale per
% species), not a taxonomy fold.
old_atlas(ismember(old_atlas.SEQ, sp_base.SEQ(sp_base.merged_SEQ == 0)), :) = [];

% Fold records for a merged_SEQ > 0 species onto their target SEQ (e.g.
% Herring Gull records -> Lesser Black-backed Gull; see the comment column on
% each such row in species_base_list.csv for whether the fold is taxonomic or
% an identification decision).
%
% Fixed 2026-09: the previous version matched old_atlas.SEQ against
% sp_base.merged_SEQ directly - but a fold-source's own SEQ (e.g. 299) is
% never itself a *value* within the merged_SEQ column (only 298, its target,
% is), so that match never found anything, and those records were then
% silently dropped once the fold-source rows were removed from sp_base below.
% 49 old-atlas squares across 6 species were affected (SEQ 299, 348, 498,
% 502, 506, 691); see git history for this file for the corrected counts.
%
% External datasets keyed to the atlas by SEQ (SABAP1-2, burns2021) need the
% same fold applied to *their* SEQ column before joining - see
% functions/fold_SEQ.m, which also knows which folds an outside dataset must
% not follow (the two that are identification decisions, not taxonomy).
fold_src = sp_base.SEQ(sp_base.merged_SEQ > 0);
fold_dst = sp_base.merged_SEQ(sp_base.merged_SEQ > 0);
[tf, loc] = ismember(old_atlas.SEQ, fold_src);
old_atlas.SEQ(tf) = fold_dst(loc(tf));

% Filter out fold-source and rejected rows from the base list: each is now
% represented under its target (or excluded entirely), not as its own concept
sp_base = sp_base(isnan(sp_base.merged_SEQ), :);

% Two namings travel side by side from here on, and nothing is overwritten:
%   common_name / scientific_name          the 1970s atlas (historical)
%   avilist_common_name / avilist_scientific_name / avilist_sort   AviList
% F_analysis.m, D_correction.m and the website's "A Bird Atlas of Kenya
% (1989)" option read the first triplet; the website's "AviList" option reads
% the second. Both are built in data/taxonomy/update_taxonomy.py, including
% the collapsed display form for a lump ("Ficedula sp."), so this file only
% has to fill the gaps that need the atlas name as a fallback.
%
% This replaces the two species that used to be renamed by hand here (the
% Fischer's Lovebird hybrid and the Ficedula group) with one general rule
% covering every multi-species concept.
%
% A concept with no resolvable AviList member falls back to the atlas name, so
% the website's AviList option always shows something. No row needs this as of
% 2026-09 (SEQ 556 was the last one, fixed by mapping it to the European +
% African Red-rumped Swallow - see data/taxonomy/SOURCES.md); the rule stays
% because a future eBird taxonomy revision can reopen the gap.
has_avilist = ~ismissing(sp_base.avilist_common_name);
sp_base.avilist_common_name(~has_avilist) = sp_base.common_name(~has_avilist);
sp_base.avilist_scientific_name(~has_avilist) = sp_base.scientific_name(~has_avilist);

% avilist_sort: AviList's own linear sequence (readtable already infers this
% numeric, NaN where AviList has no species-rank row for the concept). Falls
% back to SEQ, an imperfect but reasonable stand-in. One row needs it: SEQ 198
% African Swamphen, which AviList ranks as a subspecies while eBird splits it.
no_sort = isnan(sp_base.avilist_sort);
sp_base.avilist_sort(no_sort) = sp_base.SEQ(no_sort);

% IUCN: convert AviList's Red List code to the spelled-out category this
% pipeline's downstream figures expect (F_analysis.m hard-codes strings like
% "Critically Endangered"). A lump reports its most severe member's status
% (computed in update_taxonomy.py) so e.g. Bar-throated/Taita Apalis (SEQ
% 753, includes Critically Endangered Apalis fuscigularis) is not silently
% omitted from threat-status figures.
iucn_map = dictionary(["LC", "NT", "VU", "EN", "CR", "CR (PE)", "CR (PEW)", "EW", "EX", "DD", "NE"], ...
    ["Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered", ...
     "Critically Endangered", "Critically Endangered", "Extinct in the Wild", "Extinct", ...
     "Data Deficient", "Not Evaluated"]);
% Unset stays <missing> rather than "", so jsonencode emits null for it in
% the website export - the same "no value" spelling as birdlife_url.
sp_base.IUCN = repmat(string(missing), height(sp_base), 1);
has_iucn = ~ismissing(sp_base.iucn);
sp_base.IUCN(has_iucn) = iucn_map(sp_base.iucn(has_iucn));

% checklist_family: kept under its historical name (F_analysis.m and
% plot_tree.R both read sp.checklist_family) but now sourced from AviList
% instead of the retired 2019 checklist column of the same name.
sp_base = renamevars(sp_base, "family", "checklist_family");

% Hand-curated ecological flags, recovered 2026-09 from the last commit of
% data/species_base_list.xlsx (git history) after the CSV rewrite dropped
% them: not derivable from AviList, so update_taxonomy.py never touches
% these columns. Convert from CSV text ("true"/"false") to logical so
% jsonencode later emits real JSON booleans, not truthy strings.
for v = ["endemic", "afrotropical", "palearctic", "waterbird"]
    sp_base.(v) = sp_base.(v) == "true";
end

% Drop crosswalk-only/no-longer-needed columns: merged_SEQ (already applied
% above), avibase_id_source (internal provenance), iucn (superseded by IUCN)
sp_base = removevars(sp_base, ["merged_SEQ", "avibase_id_source", "iucn"]);

%% Format Atlas Data as Matrix

% Initialize a 3D logical array to store the presence of species in each grid cell
map_old = false(numel(g.lat), numel(g.lon), height(sp_base));

% Populate the matrix with species data
for i_sp = 1:height(sp_base)
    % Find indices of grid cells where the species is present
    id = sp_base.SEQ(i_sp) == old_atlas.SEQ;
    
    % Create a temporary grid to mark presence
    tmp = false(numel(g.lat), numel(g.lon));
    tmp(old_atlas.idg(id)) = true;
    
    % Store the species presence data in the main matrix
    map_old(:, :, i_sp) = tmp;
end

%% Visual Validation of Data

% Visualize the total number of species per grid cell
figure; 
imagesc(g.lon, g.lat, sum(map_old, 3), 'alphadata', 0.8*(sum(map_old, 3) > 0)); 
axis equal tight; 
set(gca, "YDir", "normal");
plot_google_map;  % Overlay with Google Maps
title('Number of species'); 
colorbar;

% Visualize data for a specific species
i_sp = 141;  % Example species index
figure; 
imagesc(g.lon, g.lat, map_old(:, :, i_sp), 'alphadata', 0.8*(map_old(:, :, i_sp) > 0)); 
axis equal tight; 
set(gca, "YDir", "normal");
plot_google_map;  % Overlay with Google Maps
title(sp_base.common_name(i_sp) + " (SEQ=" + sp_base.SEQ(i_sp) + ")");

%% Save Processed Data

% Save the grid data and the processed old atlas data
save('data/grid', 'g');
save('data/oldatlas', "map_old", "sp_base", "coverage_old");
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

% The digitization export preserves original SEQs (including 0a and 0b).
opts = detectImportOptions('../Digitization-of-A-Bird-Atlas-of-Kenya/export/comparison/atlas_records.csv', TextType='string');
opts = setvartype(opts, 'seq', 'string');
old_atlas = readtable('../Digitization-of-A-Bird-Atlas-of-Kenya/export/comparison/atlas_records.csv', opts);
opts_mapping = detectImportOptions('data/taxonomy/atlas_mapping.csv', TextType='string');
opts_mapping = setvartype(opts_mapping, 'atlas_seq', 'string');
mapping = readtable('data/taxonomy/atlas_mapping.csv', opts_mapping);
[found, idx] = ismember(old_atlas.seq, mapping.atlas_seq);
assert(all(found));
assert(all(old_atlas.avibase_id == mapping.atlas_avibase_id(idx)), 'Historical IDs changed: review atlas_mapping.csv.');
old_atlas.comparison_id = mapping.avibase_id(idx);
[~, old_atlas.idg] = ismember(string(old_atlas.sq_n) + old_atlas.sq_l, string(g.SqN(:)) + g.SqL(:));

% The comparison uses 1970-1984 records only (map_old). Records found only
% before 1970 or after 1984 are kept separately (map_old_other): they are not
% part of the comparison, but still show a species was once in the square
% (Kenya Bird Trends Lost & Found page).
old_atlas_other = old_atlas(ismissing(old_atlas.period_1970_1984) & old_atlas.idg > 0 & ...
    (~ismissing(old_atlas.pre_1970) | ~ismissing(old_atlas.post_1984)), :);
old_atlas = old_atlas(~ismissing(old_atlas.period_1970_1984), :);

%% Read comparison taxa and AviList names
opts_taxa = detectImportOptions('data/taxonomy/comparison_taxa.csv', TextType='string');
opts_taxa = setvartype(opts_taxa, 'atlas_seq', 'string');
sp_base = readtable('data/taxonomy/comparison_taxa.csv', opts_taxa);
assert(numel(unique(sp_base.avibase_id)) == height(sp_base));
sp_base = renamevars(sp_base, 'family', 'checklist_family');
for v = ["endemic", "afrotropical", "palearctic", "waterbird"]
    sp_base.(v) = sp_base.(v) == "true";
end
iucn_map = dictionary(["LC", "NT", "VU", "EN", "CR", "CR (PE)", "CR (PEW)", "EW", "EX", "DD", "NE"], ...
    ["Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered", ...
    "Critically Endangered", "Critically Endangered", "Extinct in the Wild", "Extinct", "Data Deficient", "Not Evaluated"]);
sp_base.IUCN = repmat(string(missing), height(sp_base), 1);
has_iucn = ~ismissing(sp_base.iucn);
sp_base.IUCN(has_iucn) = iucn_map(sp_base.iucn(has_iucn));

%% Format Atlas Data as Matrix

% Initialize a 3D logical array to store the presence of species in each grid cell
map_old = false(numel(g.lat), numel(g.lon), height(sp_base));

% Populate the matrix with species data
for i_sp = 1:height(sp_base)
    % Find indices of grid cells where the species is present
    id = sp_base.avibase_id(i_sp) == old_atlas.comparison_id;
    
    % Create a temporary grid to mark presence
    tmp = false(numel(g.lat), numel(g.lon));
    tmp(old_atlas.idg(id)) = true;
    
    % Store the species presence data in the main matrix
    map_old(:, :, i_sp) = tmp;
end

% Same for records found only outside 1970-1984
map_old_other = false(size(map_old));
for i_sp = 1:height(sp_base)
    tmp = false(numel(g.lat), numel(g.lon));
    tmp(old_atlas_other.idg(sp_base.avibase_id(i_sp) == old_atlas_other.comparison_id)) = true;
    map_old_other(:, :, i_sp) = tmp & ~map_old(:, :, i_sp);
end

%% Visual Validation of Data

% Visualize the total number of species per grid cell
figure; 
imagesc(g.lon, g.lat, sum(map_old, 3), 'alphadata', 0.8*(sum(map_old, 3) > 0)); 
axis equal tight; 
set(gca, "YDir", "normal");
if usejava('desktop'), plot_google_map; end
title('Number of species'); 
colorbar;

% Visualize data for a specific species
i_sp = 141;  % Example species index
figure; 
imagesc(g.lon, g.lat, map_old(:, :, i_sp), 'alphadata', 0.8*(map_old(:, :, i_sp) > 0)); 
axis equal tight; 
set(gca, "YDir", "normal");
if usejava('desktop'), plot_google_map; end
title(sp_base.common_name(i_sp) + " (SEQ=" + sp_base.SEQ(i_sp) + ")");

%% Save Processed Data

% Save the grid data and the processed old atlas data
save('data/grid', 'g');
save('data/oldatlas', "map_old", "map_old_other", "sp_base", "coverage_old");

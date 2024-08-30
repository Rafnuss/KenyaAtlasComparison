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

% Load the species base list
sp_base = readtable("data/species_base_list.xlsx", 'TextType', 'string');

% Remove species with MergeSEQ == 0 (i.e., species that were not merged)
old_atlas(ismember(old_atlas.SEQ, sp_base.SEQ(sp_base.merged_SEQ == 0)), :) = [];

% Merge species by replacing SEQ with merged_SEQ where applicable
[~, id] = ismember(old_atlas.SEQ, sp_base.merged_SEQ);
old_atlas.SEQ(id > 0) = sp_base.merged_SEQ(id(id > 0));

% Filter out merged species from the base list
sp_base = sp_base(isnan(sp_base.merged_SEQ), :);

% Remove unnecessary columns from the base list
sp_base = removevars(sp_base, ["merged_SEQ", "avibaseID", "inaturalistID", ...
    "observationorgID", "GBIFID", "clements_code", "SISRecID", "Min_Latitude", ...
    "Max_Latitude", "Centroid_Latitude", "Centroid_Longitude", "Range_Size"]);

% Rename specific species for consistency
id = sp_base.common_name == "Fischer's Lovebird";
sp_base.common_name(id) = "Fischer's x Yellow-collared Lovebird";
sp_base.scientific_name(id) = "Agapornis fischeri x personatus";
sp_base.IUCN(id) = missing();

% Rename Collared Flycatcher species group
id = sp_base.SEQ == 786;
sp_base.common_name(id) = "Semicollared/Pied/Collared Flycatcher";
sp_base.scientific_name(id) = "Ficedula sp.";

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
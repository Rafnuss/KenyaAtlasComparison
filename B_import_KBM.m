% Load necessary data and add paths
load('data/grid')            % Load grid data
load('data/oldatlas')        % Load old atlas data
addpath('functions/')        % Add the 'functions' directory to the MATLAB path

%% Read the curated source mapping
sp_kbm = readtable('data/kbm/species_mapping.csv', TextType='string');

%% Read GeoJSON data and convert to map matrix

% Define a finer grid for the KBM data
gn.res = 5 / 60;   % Resolution of 5 arc-minutes
gn.lon = (g.lon(1) - g.res / 2) : gn.res : (g.lon(end) + g.res / 2 - gn.res);
gn.lat = (g.lat(1) - g.res / 2) : gn.res : (g.lat(end) + g.res / 2 - gn.res);
[gn.LON, gn.LAT] = meshgrid(gn.lon, gn.lat);

% Compute pentad code for the finer grid
gn.pentad = latlon2pentad(gn.LAT, gn.LON);

% Shift grid to be cell-centered
gn.lon = gn.lon + gn.res / 2;
gn.lat = gn.lat + gn.res / 2;
[gn.LON, gn.LAT] = meshgrid(gn.lon, gn.lat);

%% Read cached source maps by KBM reference
load('data/kbm/map_kbm', 'fullp', 'adhocp', 'seq_order');
[found, idx] = ismember(sp_kbm.ref, seq_order);
assert(all(found));
fullp = fullp(:,:,idx);
adhocp = adhocp(:,:,idx);

%% Visual check of species coverage
tmp = fullp | adhocp;
tmp2 = sum(tmp, 3);
figure; hold on;
imagesc(gn.lon, gn.lat, tmp2, 'alphadata', 0.8 * (tmp2 > 0));
axis equal tight; set(gca, "YDir", "normal"); 
if usejava('desktop'), plot_google_map; end
title('Number of species')

%% Upscale map and create map_kbm
fullp(isnan(fullp)) = 0;
fullpu = blockproc(fullp, [1 1] * g.res / gn.res, @(x) sum(x.data, [1 2]));
adhocp(isnan(adhocp)) = 0;
adhocpu = blockproc(adhocp, [1 1] * g.res / gn.res, @(x) sum(x.data, [1 2]));

% Check that the coordinates are consistent (optional)
% gfu.LAT = blockproc(gn.LAT, [1 1] * g.res / gn.res, @(x) mean(x.data(:)));
% gfu.LON = blockproc(gn.LON, [1 1] * g.res / gn.res, @(x) mean(x.data(:)));
% figure; hold on;
% mesh(g.LON, g.LAT, ones(size(g.LON)), 'FaceAlpha', 0, 'EdgeColor', 'k')
% mesh(gn.LON, gn.LAT, ones(size(gn.LON)), 'FaceAlpha', 0, 'EdgeColor', 'r')
% mesh(gfu.LON, gfu.LAT, ones(size(gfu.LON)), 'FaceAlpha', 0)
% view(2); axis equal

% Compute new KBM species map
map_kbm0 = (fullpu + adhocpu) > 0;

% Visual check of the upscaled map
tmp = sum(map_kbm0, 3);
figure; hold on;
imagesc(g.lon, g.lat, tmp, 'alphadata', (tmp > 0) .* 8); 
axis equal tight; set(gca, "YDir", "normal")
if usejava('desktop'), plot_google_map; end

%% Merge maps for the same species
map_kbm = false(size(map_old));
for i_sp = 1:height(sp_base)
    id = find(sp_base.avibase_id(i_sp) == sp_kbm.avibase_id);
    map_kbm(:,:,i_sp) = any(map_kbm0(:,:,id), 3);
end

%% Get KBM effort coverage
% Check if the effort coverage data should be downloaded from the web
if false
    opts = weboptions("Timeout",60);
    
    % Download the coverage GeoJSON data from the birdmap API
    websave("data/kbm/coverage.geojson", "https://api.birdmap.africa/sabap2/v2/coverage/project/kenya?format=geoJSON")
else
    % Load and decode the pre-downloaded GeoJSON coverage data
    d = jsondecode(fileread('data/kbm/coverage.geojson'));
    
    % Initialize a matrix to store the effort coverage data
    coverage_kbm = zeros(numel(gn.lat), numel(gn.lon));
    
    % Loop through each feature (representing a pentad) in the GeoJSON file
    for i_f = 1:numel(d.features)
        % Identify the grid cells corresponding to the current pentad
        id = strcmpi(d.features(i_f).properties.pentad, gn.pentad);
        assert(sum(id(:)) > 0) % Ensure that the pentad is found in the grid
        
        % Find the indices of the grid cells
        [i2, i1] = ind2sub(size(id), find(id));
        
        % Store the full protocol total hours in the coverage matrix
        coverage_kbm(i2, i1) = d.features(i_f).properties.fullProtocol_total_hours;
        
        % Additional properties can be extracted similarly if needed:
        % fullProtocol_kbm(i2, i1) = d.features(i_f).properties.fullProtocol;
        % coverage_kbm(i2, i1, 3) = d.features(i_f).properties.adhocProtocol;
    end
    
    % Upscale the coverage matrix to match the resolution of the main grid (g)
    coverage_kbm = blockproc(coverage_kbm, [1 1] * g.res / gn.res, @(x) sum(x.data, [1 2]));
    
    % Replace NaN values with zeros in the upscaled coverage matrix
    coverage_kbm(isnan(coverage_kbm)) = 0;
end

% Use the historical atlas squares as the comparison domain.
map_kbm = map_kbm & (g.SqN > 0);
coverage_kbm(g.SqN == 0) = 0;

%% Save the final atlas data
% seq_map records which SEQ each slice of map_kbm belongs to. The three map
% layers are combined by array position downstream (map_kbm | map_ebird), so
% a layer built against a different species list silently pairs every species
% with someone else's data - which is exactly what happened when the 2026-09
% taxonomy update rebuilt oldatlas/ebirdatlas but not this file. Saving the
% key alongside the map lets D and E assert the three agree.
seq_map = sp_base.SEQ;
avibase_map = sp_base.avibase_id;
save('data/kbmatlas.mat', "map_kbm", "coverage_kbm", "seq_map", "avibase_map")

% Load necessary data and add paths
load('data/grid')            % Load grid data
load('data/oldatlas')        % Load old atlas data
addpath('functions/')        % Add the 'functions' directory to the MATLAB path

%% Get KBM species coverage

% Check if species coverage should be retrieved from the web
if false
    opts = weboptions("Timeout",60);
    
    % Download species coverage map (CSV format) from the birdmap API
    % websave("data/kbm/coverage.csv", "https://api.birdmap.africa/sabap2/v2/coverage/project/kenya/species?format=csv", opts)
    
    % Download species coverage map (GeoJSON format) from the birdmap API
    % websave("data/kbm/coverage.geojson","https://api.birdmap.africa/sabap2/v2/coverage/project/kenya?format=geoJSON")
    
    % Read species coverage data into a table
    sp_kbm = readtable("data/kbm/coverage.csv", 'TextType', 'string');
    
    % Add SEQ number to the species list by matching 'Ref' with 'ADU' in sp_base
    [~,tmp] = ismember(sp_kbm.Ref, sp_base.ADU);
    sp_kbm.SEQ(:) = nan;
    sp_kbm.SEQ(tmp > 0) = sp_base.SEQ(tmp(tmp > 0));
    
    % Save the updated species coverage data
    writetable(sp_kbm, "data/kbm/sp_kbm.csv")
    
    % Manually add missing SEQ numbers (if required)
else
    % Load pre-processed species coverage data
    sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');
end

%% Download data for all species in GeoJSON format

% Check if species data should be downloaded from the web
if false
    opts = weboptions("Timeout",60);
    
    % Loop over each species and download GeoJSON data
    for i_sp = 1:height(sp_kbm)
        filename = "data/kbm/geojson/" + sp_kbm.Ref(i_sp) + ".geojson";
        
        % Download only if species reference is not NaN
        if ~isnan(sp_kbm.Ref(i_sp))
            sp_kbm.Ref(i_sp)
            websave(filename, "https://api.birdmap.africa/sabap2/v2/summary/species/" + sp_kbm.Ref(i_sp) + "/country/kenya?format=geoJSON", opts)
            pause(5)  % Pause to avoid overloading the server
        end
    end
end


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

%% Process the species GeoJSON data
if (false)
    % Initialize matrices for storing data
    fullp = zeros(numel(gn.lat), numel(gn.lon), height(sp_kbm));
    fullpp = fullp;
    adhocp = fullp;
    
    % Loop over each species
    for i_sp = 1:height(sp_kbm)
        if ~isnan(sp_kbm.Ref(i_sp))
            % Read GeoJSON file for the species
            d = jsondecode(fileread("data/kbm/geojson/" + sp_kbm.Ref(i_sp) + ".geojson"));
            
            % Process each feature in the GeoJSON file
            for i_f = 1:numel(d.features)
                id = strcmpi(d.features(i_f).properties.pentad, gn.pentad);
                assert(sum(id(:)) > 0)
                [i2, i1] = ind2sub(size(id), find(id));
                
                % Extract data and store it in matrices
                fullp(i2, i1, i_sp) = d.features(i_f).properties.fullProtocolCards;
                fullpp(i2, i1, i_sp) = d.features(i_f).properties.fullProtocol;
                adhocp(i2, i1, i_sp) = d.features(i_f).properties.adhocProtocol;
            end
            i_sp
        end
    end
    
    % Save processed data
    seq_order = sp_kbm.Ref;
    save('data/kbm/map_kbm', "fullp", "fullpp", "adhocp", "seq_order", "gn")
else
    % Load pre-processed KBM map data
    load('data/kbm/map_kbm')
    assert(all(seq_order == sp_kbm.Ref))
end

%% Visual check of species coverage
tmp = fullp | adhocp;
tmp2 = sum(tmp, 3);
figure; hold on;
imagesc(gn.lon, gn.lat, tmp2, 'alphadata', 0.8 * (tmp2 > 0));
axis equal tight; set(gca, "YDir", "normal"); 
plot_google_map;
title('Number of species')

% Example: Visualize coverage for a specific species
i_sp = find(strcmp(sp_kbm.Common_species, "Red-footed"));
tmp2 = sum(tmp(:,:,i_sp), 3);
figure; hold on;
imagesc(gn.lon, gn.lat, tmp2, 'alphadata', 0.8 * (tmp2 > 0)); 
axis equal tight; set(gca, "YDir", "normal"); 
plot_google_map;
title(sp_kbm.Common_species(i_sp) + " " + sp_kbm.Common_group(i_sp) + " (ADU:" + sp_kbm.Ref(i_sp) + " | SEQ:" + sp_kbm.SEQ(i_sp) + ")")

% Compute number of pentads recorded for each species (if needed)
% sp_kbm.nb_pentad_recorded = squeeze(sum(tmp, [1 2]));

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
plot_google_map;

%% Merge maps for the same species
map_kbm = false(size(map_old));
for i_sp = 1:height(sp_base)
    id = find(sp_base.SEQ(i_sp) == sp_kbm.SEQ);
    map_kbm(:,:,i_sp) = any(map_kbm0(:,:,id), 3);
end

% Identify species not kept in the final map
disp("Species not kept: Should all be entry error in KBM.")
sp_kbm(~ismember(sp_kbm.SEQ, sp_base.SEQ) & sp_kbm.Number_pentad_recorded > 0, :)

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

%% Save the final atlas data
% Save the KBM species map and effort coverage data to a MAT-file
save('data/kbmatlas.mat', "map_kbm", "coverage_kbm")

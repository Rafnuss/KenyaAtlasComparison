% Load necessary grid and atlas data
load('data/grid')
load('data/oldatlas')

%% Read and filter eBird data

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

%% Compute coverage map

% Summarize eBird data by checklist
ebd_checklist = groupsummary(ebd, {'SAMPLINGEVENTIDENTIFIER', 'ALLSPECIESREPORTED', 'DURATIONMINUTES', 'EFFORTDISTANCEKM', 'PROTOCOLTYPE', 'NUMBEROBSERVERS', 'idg'});

% Summarize data by grid cell
ebd_grid = groupsummary(ebd_checklist, "idg", "sum", "DURATIONMINUTES");

% Initialize coverage map and fill with data
coverage_ebird = nan(size(g.LAT));
coverage_ebird(ebd_grid.idg) = ebd_grid.sum_DURATIONMINUTES / 60; % Convert minutes to hours
coverage_ebird(isnan(coverage_ebird)) = 0;

%% Combine species/grid data

% Select and simplify eBird data for further processing
ebd = table(ebd.LATITUDE, ebd.LONGITUDE, ebd.COMMONNAME, ebd.SCIENTIFICNAME, ebd.CATEGORY, ebd.idg, 'VariableNames', {'lat', 'lon', 'common_name', 'scientific_name', 'category', 'idg'});
ebd = unique(ebd, "sorted");

%% Taxonomy Matching

% Load the species taxonomy data for eBird
sp_ebird = readtable('data/eBird/sp_ebird.xlsx', 'TextType', 'string');

% Match species in eBird data with those in the existing species base
[Lia, Locb] = ismember(ebd.scientific_name, sp_ebird.scientific_name);

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

% Save the processed eBird map and coverage data to a MAT-file
save('data/ebirdatlas.mat', "map_ebird", "coverage_ebird")
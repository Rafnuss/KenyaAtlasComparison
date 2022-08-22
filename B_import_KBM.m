
load('data/grid')
load('data/oldatlas')
addpath('functions/')

%% Get kbm coverage
if false
    opts = weboptions("Timeout",60);
    
    % get coverage map
    % websave("data/kbm/coverage.csv","https://api.birdmap.africa/sabap2/v2/coverage/project/kenya/species?format=csv", opts)
    sp_kbm = readtable("data/kbm/coverage.csv", 'TextType', 'string');
    
    % Add SEQ number to the species list
    [~,tmp] = ismember(sp_kbm.Ref,sp_old.ADU); 
    sp_kbm.SEQ(:) = nan;
    sp_kbm.SEQ(tmp>0) = sp_old.SEQ(tmp(tmp>0));
    writetable(sp_kbm,"data/kbm/sp_kbm.csv")

    % Manuall add missing SEQ number
else
    sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');
end


%% Download the data for all species in geojson
if false
    opts = weboptions("Timeout",60);
    for i_sp=1:height(sp_kbm)
        filename = "data/kbm/geojson/"+sp_kbm.Ref(i_sp)+".geojson";
        if ~isnan(sp_kbm.Ref(i_sp)) && ~exist(filename,'file')
            sp_kbm.Ref(i_sp)
            websave(filename,"https://api.birdmap.africa/sabap2/v2/summary/species/"+sp_kbm.Ref(i_sp)+"/country/kenya?format=geoJSON",opts)
            pause(5)
        end
    end
end

%% Read geojson and convert to maptrix

% define grid
gn.res = 5 / 60;
gn.lon = (g.lon(1)-g.res/2):gn.res:(g.lon(end)+g.res/2-gn.res);
gn.lat = (g.lat(1)-g.res/2):gn.res:(g.lat(end)+g.res/2-gn.res);
[gn.LAT,gn.LON] = meshgrid(gn.lat,gn.lon);

% compute pentad code
gn.pentad = latlon2pentad(gn.LAT, gn.LON);

% cell-centered
gn.lon = gn.lon + gn.res/2;
gn.lat = gn.lat + gn.res/2;
[gn.LAT,gn.LON] = meshgrid(gn.lat,gn.lon);

%% read json data
if (false)
    fullp = zeros(numel(gn.lat), numel(gn.lon), height(sp_kbm));
    adhocp = zeros(numel(gn.lat), numel(gn.lon), height(sp_kbm));
    
    for i_sp=1:height(sp_kbm)
        if ~isnan(sp_kbm.Ref(i_sp))
            d = loadjson("data/kbm/geojson/"+sp_kbm.Ref(i_sp)+".geojson");
        
            for i_f = 1:numel(d.features)
                id = strcmpi(d.features{i_f}.properties.pentad, gn.pentad);
                assert(sum(id(:))>0)
                [i1,i2]=ind2sub(size(id),find(id));
                fullp(i2,i1,i_sp) = d.features{i_f}.properties.full_0x20_protocol_0x20_cards;
                adhocp(i2,i1,i_sp) = d.features{i_f}.properties.adhoc_0x20_protocol;
            end
            i_sp
        end
    end
    % Save new atlas
    save('data/kbm/map_kbm',"fullp","adhocp")
end

%% 
load('data/kbm/map_kbm')

% Visual check
tmp = fullp | adhocp;
tmp2 = sum(tmp,3);
figure; hold on;
imagesc(gn.lon,gn.lat,tmp2,'alphadata',0.8*(tmp2>0)); 
axis equal tight; set(gca,"YDir","normal"); plot_google_map;
title('Number of species')

i_sp = 140;
tmp2 = sum(tmp(:,:,i_sp),3);
figure; hold on;
imagesc(gn.lon,gn.lat,tmp2,'alphadata',0.8*(tmp2>0)); 
axis equal tight; set(gca,"YDir","normal"); plot_google_map;
title(sp_kbm.Common_species(i_sp)+" " + sp_kbm.Common_group(i_sp)+" ("+sp_kbm.Ref(i_sp)+")")



%% Upscale map and create map_kbm
fullp(isnan(fullp))=0;
fullpu = blockproc(fullp,[1 1]*g.res/gn.res,@(x) sum(x.data,[1 2]));
adhocp(isnan(adhocp))=0;
adhocpu = blockproc(adhocp,[1 1]*g.res/gn.res,@(x) sum(x.data,[1 2]));

% check coordinate (should be the same as g)
%gfu.LAT = blockproc(gn.LAT,[1 1]*g.res/gn.res,@(x) mean(x.data(:)));
%gfu.LON = blockproc(gn.LON,[1 1]*g.res/gn.res,@(x) mean(x.data(:)));
% figure; hold on;
% mesh(g.LON,g.LAT,ones(size(g.LON)),'FaceAlpha',0,'EdgeColor','k')
% mesh(gn.LON,gn.LAT,ones(size(gn.LON)),'FaceAlpha',0,'EdgeColor','r')
% mesh(gfu.LON,gfu.LAT,ones(size(gfu.LON)),'FaceAlpha',0)
% view(2); axis equal

% compute new map
map_kbm0 = (fullpu+adhocpu) > 0;

% Visual check
tmp = sum(map_kbm0,3);
figure; hold on;
imagesc(g.lon, g.lat, tmp,'alphadata',(tmp>0).*8); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map;


%% Merge maps for the same species
map_kbm = false(size(map_old));
for i_sp=1:height(sp_old)
    id = find(sp_old.SEQ(i_sp)==sp_kbm.SEQ);
    map_kbm(:,:,i_sp) = any(map_kbm0(:,:,id),3);
end


%% Save
save('data/kbmatlas.mat',"map_kbm")
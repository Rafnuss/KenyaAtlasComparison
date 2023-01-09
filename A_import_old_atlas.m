
%% 1. Get grid with SqN and SqL

% Read grid data generated with the R code
grid = jsondecode(fileread('data/oldatlas/grid.geojson'));

% retrieve the coordinate of the center of each cell
coord=nan(2,numel(grid.features));
for i_f =1:numel(grid.features)
    coord(:,i_f) = round(squeeze(mean(grid.features(i_f).geometry.coordinates)),1);
end

% Define the grid
g.res = min(diff(sort(unique(coord(1,:)))));
g.lon =  round(min(coord(1,:)):g.res:max(coord(1,:)),1);
g.lat =  round(min(coord(2,:)):g.res:max(coord(2,:)),1);
[g.LON,g.LAT] = meshgrid(g.lon,g.lat);

g.SqL = strings(numel(g.lat), numel(g.lon));
g.SqN = zeros(numel(g.lat), numel(g.lon));
coverage_old = strings(numel(g.lat), numel(g.lon));
for i_f = 1:numel(grid.features)
    id_lat = coord(2,i_f)==g.lat;
    id_lon = g.lon==coord(1,i_f);
    assert(sum(id_lat)==1 & sum(id_lon)==1)
    g.SqL(id_lat,id_lon) = grid.features(i_f).properties.SqL;
    g.SqN(id_lat,id_lon) = grid.features(i_f).properties.SqN;
    coverage_old(id_lat,id_lon) = grid.features(i_f).properties.coverage;
end


cov = nan(size(coverage_old));
cov(coverage_old=="0")=0;
cov(coverage_old=="1-10")=1;
cov(coverage_old=="11-30")=2;
cov(coverage_old=="31-50")=3;
cov(coverage_old=="51-75")=4;
cov(coverage_old=="75-100")=5;
figure; imagesc(cov); set(gca,'ydir','normal')

%% Import atlas data
old_atlas = readtable("data/oldatlas/A Bird Atlas of Kenya_v5.xlsx", 'TextType', 'string');

% Keep only the data of interest
old_atlas = old_atlas(:,["SEQ" "SqN" "SqL" "pre_1970" "x1970_1984"]);

% Remove recent data (post 1984)
old_atlas(ismissing(old_atlas.pre_1970) & ismissing(old_atlas.x1970_1984),:)=[];

% Match the grid
[~,old_atlas.idg]=ismember(string(old_atlas.SqN)+old_atlas.SqL,string(g.SqN(:))+g.SqL(:));



%% Script to create the dataset to compare the old and new Kenyan bird atlas

sp_old = readtable("data/oldatlas/A Bird Atlas of Kenya_base_list.xlsx", 'TextType', 'string');

% Delete species with MergeSEQ==0
sp_old(sp_old.MergedSEQ==0,:)
old_atlas(ismember(old_atlas.SEQ,sp_old.SEQ(sp_old.MergedSEQ==0)),:)=[];

% Merge species by replacing SEQ for MergedSEQ>0
[~,id] = ismember(old_atlas.SEQ, sp_old.MergedSEQ);
old_atlas.SEQ(id>0) = sp_old.MergedSEQ(id(id>0));

sp_old = sp_old(isnan(sp_old.MergedSEQ),:);


%% Format as matrix
map_old = false(numel(g.lat), numel(g.lon), height(sp_old));

for i_sp = 1:height(sp_old)
    id = sp_old.SEQ(i_sp) == old_atlas.SEQ;
    tmp = false(numel(g.lat), numel(g.lon));
    tmp(old_atlas.idg(id))=true;
    map_old(:,:,i_sp) = tmp;
end

%% Check visually that everything is correct.
figure; imagesc(g.lon,g.lat,sum(map_old,3),'alphadata',0.8*(sum(map_old,3)>0)); axis equal tight; set(gca,"YDir","normal")
plot_google_map; title('Number of species'); colorbar;

i_sp = 141;
figure; imagesc(g.lon,g.lat,map_old(:,:,i_sp),'alphadata',0.8*(map_old(:,:,i_sp)>0)); axis equal tight; set(gca,"YDir","normal")
plot_google_map; title(sp_old.CommonName(i_sp)+" (SEQ="+sp_old.SEQ(i_sp)+")")



%% Save
save('data/grid','g')
save('data/oldatlas',"map_old","sp_old","coverage_old")
%% Import data

addpath('functions/')
load('data/grid_corr')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
% grid = loadjson('data/oldatlas/grid.geojson');
grid = jsondecode(fileread("data/oldatlas/grid.geojson"));

%%

map_new = map_kbm | map_ebird;
map_data = grid;
for i_g=1:numel(map_data.features)

    % Change polygon to marker
    map_data.features(i_g).geometry.type = "Point";
    coord = squeeze(map_data.features(i_g).geometry.coordinates);
    map_data.features(i_g).geometry.coordinates = fliplr(mean(coord(1:4,:)));

    id = map_data.features(i_g).properties.SqL==g.SqL & map_data.features(i_g).properties.SqN==g.SqN;
    assert(sum(id(:))==1)
    id_old = map_old(repmat(id,1,1,size(map_old,3)));
    id_new = map_new(repmat(id,1,1,size(map_old,3)));
    prop = map_data.features(i_g).properties;
    prop.nb_lkgd = [sum(id_old & ~id_new) sum(id_old & id_new) sum(~id_old & id_new) sum(id_new-id_old)];
    if sum(id_old)>1
        prop.SEQ_old = sp_base.SEQ(id_old);
    else
        prop.SEQ_old = {sp_base.SEQ(id_old)};
    end
    if sum(id_new)>1
        prop.SEQ_new = sp_base.SEQ(id_new);
    else
        prop.SEQ_new = {sp_base.SEQ(id_new)};
    end
    prop.Sq = string(prop.SqN) + prop.SqL;

    prop = rmfield(prop,{'SqN','SqL', 'coverage'});

    prop.cov_new = round(g.coverage_new(id));
    prop.cov_old = g.coverage_old(id);
    prop.mask =  g.mask(id);
    prop.corr = g.corr(id);
    map_data.features(i_g).properties = prop;
end

fid = fopen('export/website/map_data.json','w');
fprintf(fid,'%s',jsonencode(map_data));
fclose(fid);

%% Export grid
grid_web = grid;
for i_g=1:numel(grid_web.features)
    grid_web.features(i_g).properties = struct();
    grid_web.features(i_g).geometry.coordinates = grid_web.features(i_g).geometry.coordinates;
end

fid = fopen('export/website/grid.json','w');
fprintf(fid,'%s',jsonencode(grid_web));
fclose(fid);

%% SP_old
sp_base2=sp_base;
% sp_base2.IUCN(sp_base2.IUCN=="Critically Endangered")="CR"  ;
% sp_base2.IUCN(sp_base2.IUCN=="Endangered")="EN"  ;
% sp_base2.IUCN(sp_base2.IUCN=="Vulnerable")="VU"  ;
% sp_base2.IUCN(sp_base2.IUCN=="Near Threatened")="NT"  ;
% sp_base2.IUCN(sp_base2.IUCN=="Least Concern")="LC"  ;
% sp_base2.IUCN(sp_base2.IUCN=="Data Deficient")="DD"  ;


sp_base2.nb_lkgd = [
    reshape(sum(map_old & ~map_new,[1 2]),[],1)... %lost
    reshape(sum(map_old & map_new,[1 2]),[],1)... %kept
    reshape(sum(~map_old & map_new,[1 2]),[],1) ... %gain
    reshape(sum(map_new-map_old,[1 2]),[],1)]; %diff

sp_base2.per_lkgd = sp_base2.nb_lkgd ./ sum(sp_base2.nb_lkgd(:,1:3),2);

mask3d = ~repmat(g.mask, 1, 1, size(map_old,3));
sp_base2.nb_lkgd_gc = [
    reshape(sum(map_old & ~map_new & mask3d,[1 2]),[],1)...
    reshape(sum(map_old & map_new & mask3d,[1 2]),[],1)... 
    reshape(sum(~map_old & map_new & mask3d,[1 2]),[],1) ...
    reshape(sum(map_new-map_old & mask3d,[1 2]),[],1)];

sp_base2.per_lkgd_gc = sp_base2.nb_lkgd_gc ./ sum(sp_base2.nb_lkgd_gc(:,1:3),2);
sp_base2.per_lkgd_gc(isnan(sp_base2.per_lkgd_gc)) = 0;

sp_base2 = sortrows(sp_base2,"SEQ");
sp_base2 = removevars(sp_base2,["comment","checklist_family", "mass", "habitat", "ADU","migration",...
    "habitat_density", "trophic_level", "trophic_niche", "primary_lifestyle"]);

sp_ebird = readtable('data/eBird/sp_ebird.xlsx');
sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');

sp_base2.kbm(:)=cell(1);
sp_base2.ebird(:)=cell(1);
for i_sp=1:height(sp_base2)
    % KBM
    tmp = sp_kbm.Ref(sp_base2.SEQ(i_sp) ==sp_kbm.SEQ);
    if numel(tmp)==1, tmp = {tmp}; end
    sp_base2.kbm{i_sp} = tmp;
    % eBird
    tmp = sp_ebird.species_code(sp_base2.SEQ(i_sp) ==sp_ebird.SEQ);
    if numel(tmp)==1, tmp = {tmp}; end
    sp_base2.ebird{i_sp} = tmp;
end

fid = fopen('export/website/sp_base.json','w');
fprintf(fid,'%s',jsonencode(sp_base2));
fclose(fid);

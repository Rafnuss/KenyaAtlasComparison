%% Import data

load('data/grid')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
grid = loadjson('data/oldatlas/grid.geojson');

%% Export target data

map_new = map_kbm | map_ebird;
grid_target = grid;
for i_g=1:numel(grid_target.features)
     id = grid_target.features{i_g}.properties.SqL==g.SqL & grid_target.features{i_g}.properties.SqN==g.SqN;
     assert(sum(id(:))==1)
     id_old = map_old(repmat(id,1,1,size(map_old,3)));
     id_new = map_new(repmat(id,1,1,size(map_old,3)));
     prop = grid_target.features{i_g}.properties;
     prop.nb_lkgd = [sum(id_old & ~id_new) sum(id_old & id_new) sum(~id_old & id_new) sum(id_new-id_old)];
     if sum(id_old)>1
        prop.SEQ_old = sp_old.SEQ(id_old);
     else
        prop.SEQ_old = {sp_old.SEQ(id_old)};
     end
     if sum(id_new)>1
        prop.SEQ_new = sp_old.SEQ(id_new);
     else
        prop.SEQ_new = {sp_old.SEQ(id_new)};
     end
     prop.Sq = string(prop.SqN) + prop.SqL;
     prop = rmfield(prop,{'SqN','SqL','coverage'});
     grid_target.features{i_g}.properties = prop;
     grid_target.features{i_g}.geometry.coordinates = {grid_target.features{i_g}.geometry.coordinates};
end

fid = fopen('export/website/grid_target.json','w');
fprintf(fid,'%s',jsonencode(grid_target));
fclose(fid);


%% SP_old
sp_old.nb_lkgd = [
    reshape(sum(map_old & ~map_new,[1 2]),[],1)...
    reshape(sum(map_old & map_new,[1 2]),[],1)... 
    reshape(sum(~map_old & map_new,[1 2]),[],1) ...
    reshape(sum(map_new-map_old,[1 2]),[],1)];
% sp_old = rmfield(sp_old,{'MergedSEQ','sort_2019'});

sp_old=sortrows(sp_old,"SEQ");

sp_ebird = readtable('data/eBird/sp_ebird.xlsx');

sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');

sp_old.kbm(:)="";
sp_old.ebird(:)="";
for i_sp=1:height(sp_old)
    % KBM
    sp_old.kbm(i_sp) = join(string(sp_kbm.Ref(sp_old.SEQ(i_sp) ==sp_kbm.SEQ)),',');
    % eBird
    sp_old.ebird(i_sp) = join(string(sp_ebird.species_code(sp_old.SEQ(i_sp) ==sp_ebird.SEQ)),',');
end

fid = fopen('export/website/sp_old.json','w');
fprintf(fid,'%s',jsonencode(sp_old));
fclose(fid);

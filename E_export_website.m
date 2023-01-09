%% Import data

load('data/grid_corr')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
grid = loadjson('data/oldatlas/grid.geojson');

%%

map_new = map_kbm | map_ebird;
map_data = grid;
for i_g=1:numel(map_data.features)

    % Change polygon to marker
    map_data.features{i_g}.geometry.type = "Point";
    map_data.features{i_g}.geometry.coordinates = fliplr(mean(map_data.features{i_g}.geometry.coordinates(1:4,:)));

    id = map_data.features{i_g}.properties.SqL==g.SqL & map_data.features{i_g}.properties.SqN==g.SqN;
    assert(sum(id(:))==1)
    id_old = map_old(repmat(id,1,1,size(map_old,3)));
    id_new = map_new(repmat(id,1,1,size(map_old,3)));
    prop = map_data.features{i_g}.properties;
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

    prop = rmfield(prop,{'SqN','SqL', 'coverage'});

    prop.coverage_new = g.coverage_new(id);
    prop.coverage_old = g.coverage_old(id);
    prop.mask = g.mask(id);
    prop.corr = g.corr(id);
    map_data.features{i_g}.properties = prop;
end

fid = fopen('export/website/map_data.json','w');
fprintf(fid,'%s',jsonencode(map_data));
fclose(fid);



%% SP_old
sp_old2=sp_old;
sp_old2.IUCN(sp_old2.IUCN=="Critically Endangered")="CR"  ;
sp_old2.IUCN(sp_old2.IUCN=="Endangered")="EN"  ;
sp_old2.IUCN(sp_old2.IUCN=="Vulnerable")="VU"  ;
sp_old2.IUCN(sp_old2.IUCN=="Near Threatened")="NT"  ;
sp_old2.IUCN(sp_old2.IUCN=="Least Concern")="LC"  ;
sp_old2.IUCN(sp_old2.IUCN=="Data Deficient")="DD"  ;


sp_old2.nb_lkgd = [
    reshape(sum(map_old & ~map_new,[1 2]),[],1)...
    reshape(sum(map_old & map_new,[1 2]),[],1)... 
    reshape(sum(~map_old & map_new,[1 2]),[],1) ...
    reshape(sum(map_new-map_old,[1 2]),[],1)];

sp_old2.per_lkgd = sp_old2.nb_lkgd ./ sum(sp_old2.nb_lkgd,2);

sp_old2 = sortrows(sp_old2,"SEQ");
sp_old2 = sp_old2(~isnan(sp_old2.SEQ) & isnan(sp_old2.MergedSEQ),:);
sp_old2 = removevars(sp_old2,{'MergedSEQ','sort_2019'});

sp_ebird = readtable('data/eBird/sp_ebird.xlsx');

sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');

sp_old2.kbm(:)=cell(1);
sp_old2.ebird(:)=cell(1);
for i_sp=1:height(sp_old2)
    % KBM
    tmp = sp_kbm.Ref(sp_old2.SEQ(i_sp) ==sp_kbm.SEQ);
    if numel(tmp)==1, tmp = {tmp}; end
    sp_old2.kbm{i_sp} = tmp;
    % eBird
    tmp = sp_ebird.species_code(sp_old2.SEQ(i_sp) ==sp_ebird.SEQ);
    if numel(tmp)==1, tmp = {tmp}; end
    sp_old2.ebird{i_sp} = tmp;
end

fid = fopen('export/website/sp_old.json','w');
fprintf(fid,'%s',jsonencode(sp_old2));
fclose(fid);



%% Export target data : NOT USED ANYMORE

map_new = map_kbm | map_ebird;
map_data = grid;
for i_g=1:numel(map_data.features)
     id = map_data.features{i_g}.properties.SqL==g.SqL & map_data.features{i_g}.properties.SqN==g.SqN;
     assert(sum(id(:))==1)
     id_old = map_old(repmat(id,1,1,size(map_old,3)));
     id_new = map_new(repmat(id,1,1,size(map_old,3)));
     prop = map_data.features{i_g}.properties;
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
     map_data.features{i_g}.properties = prop;
     map_data.features{i_g}.geometry.coordinates = {map_data.features{i_g}.geometry.coordinates};
end

fid = fopen('export/website/grid_target.json','w');
fprintf(fid,'%s',jsonencode(map_data));
fclose(fid);

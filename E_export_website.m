%% Import data

addpath('functions/')
load('data/grid_corr')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
% grid = loadjson('data/oldatlas/grid.geojson');
grid = jsondecode(fileread("data/oldatlas/grid.geojson"));

% map_new below is map_kbm | map_ebird, combined by array position, so all
% three layers must share one species list - see
% functions/check_seq_alignment.m for what goes wrong when they don't.
check_seq_alignment(sp_base, "data/kbmatlas.mat", "data/ebirdatlas.mat")

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
    id_other = map_old_other(repmat(id,1,1,size(map_old,3)));
    prop = map_data.features(i_g).properties;
    prop.nb_lkgd = [sum(id_old & ~id_new) sum(id_old & id_new) sum(~id_old & id_new) sum(id_new-id_old)];
    % Cell arrays keep zero-, one- and many-taxon memberships flat in JSON.
    prop.taxa_old = cellstr(sp_base.avibase_id(id_old));
    prop.taxa_new = cellstr(sp_base.avibase_id(id_new));
    % Atlas records only before 1970 or after 1984: not compared, but historical.
    prop.taxa_old_other = cellstr(sp_base.avibase_id(id_other));
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

sp_base2.nb_lkgd = [
    reshape(sum(map_old & ~map_new,[1 2]),[],1)... %lost
    reshape(sum(map_old & map_new,[1 2]),[],1)... %kept
    reshape(sum(~map_old & map_new,[1 2]),[],1) ... %gain
    reshape(sum(map_new-map_old,[1 2]),[],1)]; %diff

sp_base2.per_lkgd = sp_base2.nb_lkgd ./ sum(sp_base2.nb_lkgd(:,1:3),2);
sp_base2.per_lkgd(isnan(sp_base2.per_lkgd)) = 0;

mask3d = ~repmat(g.mask, 1, 1, size(map_old,3));
sp_base2.nb_lkgd_gc = [
    reshape(sum(map_old & ~map_new & mask3d,[1 2]),[],1)...
    reshape(sum(map_old & map_new & mask3d,[1 2]),[],1)... 
    reshape(sum(~map_old & map_new & mask3d,[1 2]),[],1) ...
    reshape(sum((map_new-map_old).*mask3d,[1 2]),[],1)];

sp_base2.per_lkgd_gc = sp_base2.nb_lkgd_gc ./ sum(sp_base2.nb_lkgd_gc(:,1:3),2);
sp_base2.per_lkgd_gc(isnan(sp_base2.per_lkgd_gc)) = 0;

sp_base2 = sortrows(sp_base2,"SEQ");

% Website naming and source links come from the same explicit mappings.
sp_base2 = renamevars(sp_base2, "checklist_family", "family");
sp_base2 = sp_base2(:, ["SEQ", "avibase_id", "common_name", "scientific_name", ...
    "avilist_common_name", "avilist_scientific_name", "avilist_sort", "avilist_members", "name_source", "atlas_seq", ...
    "family", "order", "IUCN", "birdlife_url", "endemic", "afrotropical", "palearctic", "waterbird", ...
    "flag", "nb_lkgd", "per_lkgd", "nb_lkgd_gc", "per_lkgd_gc"]);
sp_kbm = readtable('data/kbm/species_mapping.csv', TextType='string');
sp_ebird = readtable('data/eBird/species_mapping.csv', TextType='string');
atlas_options = detectImportOptions('data/taxonomy/atlas_mapping.csv', TextType='string');
atlas_options = setvartype(atlas_options, 'atlas_seq', 'string');
atlas_mapping = readtable('data/taxonomy/atlas_mapping.csv', atlas_options);
taxon_options = detectImportOptions('../Digitization-of-A-Bird-Atlas-of-Kenya/data/curated/atlas_taxonomy.csv', TextType='string');
taxon_options = setvartype(taxon_options, 'SEQ', 'string');
atlas_taxa = readtable('../Digitization-of-A-Bird-Atlas-of-Kenya/data/curated/atlas_taxonomy.csv', taxon_options);
sp_base2.atlas_concepts = cell(height(sp_base2), 1);
sp_base2.kbm = cell(height(sp_base2), 1);
sp_base2.ebird = cell(height(sp_base2), 1);
for i_sp = 1:height(sp_base2)
    id = sp_base2.avibase_id(i_sp);
    atlas_rows = atlas_mapping(atlas_mapping.avibase_id == id,:);
    [~, atlas_idx] = ismember(string(atlas_rows.atlas_seq),string(atlas_taxa.SEQ));
    sp_base2.atlas_concepts{i_sp} = num2cell(table2struct(table(atlas_rows.atlas_avibase_id, atlas_taxa.common_name(atlas_idx), atlas_taxa.scientific_name(atlas_idx), 'VariableNames', {'avibase_id','common_name','scientific_name'})));
    sp_base2.kbm{i_sp} = num2cell(unique(sp_kbm.ref(sp_kbm.avibase_id == id)));
    codes = sp_ebird.link_code(sp_ebird.avibase_id == id);
    sp_base2.ebird{i_sp} = cellstr(unique(codes(~ismissing(codes))));
end
fid = fopen('export/website/sp_base.json','w');
fprintf(fid,'%s',jsonencode(sp_base2));
fclose(fid);

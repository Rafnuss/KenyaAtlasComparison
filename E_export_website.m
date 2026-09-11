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
    % num2cell so jsonencode always emits a flat array, whatever the count.
    % The previous branch wrapped in a plain cell, which is right for one
    % species ({5} -> [5]) but produced [[]] for a square with none - and the
    % website flattens one level, so that nested [] landed in the SEQ set as
    % an array object: it never matches a SEQ, inflating the square's count
    % and dropping a species. 66 squares emitted a nested SEQ_new (2026-09).
    prop.SEQ_old = num2cell(sp_base.SEQ(id_old));
    prop.SEQ_new = num2cell(sp_base.SEQ(id_new));
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

mask3d = ~repmat(g.mask, 1, 1, size(map_old,3));
sp_base2.nb_lkgd_gc = [
    reshape(sum(map_old & ~map_new & mask3d,[1 2]),[],1)...
    reshape(sum(map_old & map_new & mask3d,[1 2]),[],1)... 
    reshape(sum(~map_old & map_new & mask3d,[1 2]),[],1) ...
    reshape(sum(map_new-map_old & mask3d,[1 2]),[],1)];

sp_base2.per_lkgd_gc = sp_base2.nb_lkgd_gc ./ sum(sp_base2.nb_lkgd_gc(:,1:3),2);
sp_base2.per_lkgd_gc(isnan(sp_base2.per_lkgd_gc)) = 0;

sp_base2 = sortrows(sp_base2,"SEQ");

% The exported JSON is the contract the website reads (asserted by
% tests/dataSchema.test.js in Rafnuss/KenyaBirdTrends), so emit exactly the
% fields it needs under their public names:
%   - checklist_family is an internal alias kept only so F_analysis.m and
%     plot_tree.R keep working; the website sees it as `family`
%   - ebird_code/avilist_members are the "|"-joined machine-readable forms,
%     superseded here by the `ebird` array and `avilist_scientific_name`
%   - ADU (SAFRING code) and comment are internal-only
sp_base2 = renamevars(sp_base2, "checklist_family", "family");
sp_base2 = removevars(sp_base2, ["comment", "ADU"]);

% AVONET traits stay out of the export: the site filters on the four
% ecological flags but never reads these, and they are ~15% of the payload.
% F_analysis.m is what uses them (see data/taxonomy/SOURCES.md).
sp_base2 = removevars(sp_base2, ["mass", "habitat", "habitat_density", ...
    "migration", "trophic_level", "trophic_niche", "primary_lifestyle", "range_size"]);

sp_kbm = readtable("data/kbm/sp_kbm.xlsx", 'TextType', 'string');

sp_base2.kbm(:)=cell(1);
sp_base2.ebird(:)=cell(1);
for i_sp=1:height(sp_base2)
    % KBM
    tmp = sp_kbm.Ref(sp_base2.SEQ(i_sp) ==sp_kbm.SEQ);
    if numel(tmp)==1, tmp = {tmp}; end
    sp_base2.kbm{i_sp} = tmp;
end

% eBird codes: split from the crosswalk's ebird_code (data/taxonomy/
% update_taxonomy.py), not re-queried from sp_ebird.xlsx by SEQ as before -
% that pulled in every eBird taxon mapped to this SEQ, including slash/spuh
% entries (e.g. SEQ 786 got 5 codes: y00820, ficedu1, colfly1, eupfly1,
% semfly1). Only the last 3 are real species with an eBird species page;
% the other 2 produced dead "eBird-y00820" links on the site. ebird_code
% already carries only the species-level codes.
no_code = ismissing(sp_base2.ebird_code) | sp_base2.ebird_code == "";
sp_base2.ebird_code(no_code) = "";
sp_base2.ebird = arrayfun(@(s) cellstr(split(s, "|")), sp_base2.ebird_code, 'UniformOutput', false);
% SEQ 556 (Red-rumped Swallow) has no resolved eBird code (see SOURCES.md):
% split("") gives {''} rather than {}, which would render one blank button.
% Note readtable gives <missing>, not "", for a blank CSV cell - test both.
sp_base2.ebird(no_code) = {cell(1, 0)};

% Now that `ebird` carries the codes, drop the "|"-joined machine-readable
% forms: the website contract uses the array and avilist_scientific_name.
sp_base2 = removevars(sp_base2, ["ebird_code", "avilist_members"]);

fid = fopen('export/website/sp_base.json','w');
fprintf(fid,'%s',jsonencode(sp_base2));
fclose(fid);

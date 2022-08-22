
load('data/grid')
load('data/oldatlas')

%% Read and filter ebird data

ebd = readtable("data/eBird/ebd_KE_prv_relMar-2021/ebd_KE_prv_relMar-2021.txt");

ebd = ebd(ebd.CATEGORY=="species" & year(ebd.OBSERVATIONDATE)>2000,:); 

[~,id_lat]=min((g.lat-ebd.LATITUDE).^2,[],2);
[~,id_lon]=min((g.lon-ebd.LONGITUDE).^2,[],2);
ebd.idg = sub2ind(size(g.LAT),id_lon,id_lat);

%% explore
if false
    groupsummary(ebd,"PROTOCOLTYPE")
    
    numel(unique(ebd.COMMONNAME))
    
    tmp = groupsummary(ebd,"OBSERVERID");
    height(tmp)
    sum(tmp.GroupCount>100)
    
    
    ebd_checklist = groupsummary(ebd,{'ALLSPECIESREPORTED','DURATIONMINUTES','EFFORTDISTANCEKM','LATITUDE','LONGITUDE','PROTOCOLTYPE','NUMBEROBSERVERS','id_grid'});
    
    sum(ebd_checklist.ALLSPECIESREPORTED)
    
    
    tmp2 = groupsummary(ebd_checklist(ebd_checklist.ALLSPECIESREPORTED==1,:),"id_grid","sum",["DURATIONMINUTES",'EFFORTDISTANCEKM','NUMBEROBSERVERS']);
    tmp=nan(size(g.LAT));
    tmp(tmp2.idg) = tmp2.sum_DURATIONMINUTES;
    
    figure; imagesc(g.lon,g.lat,tmp)
end
%%
ebd = table(ebd.LATITUDE,ebd.LONGITUDE,ebd.COMMONNAME,ebd.SCIENTIFICNAME, ebd.idg,'VariableNames',{'lat','lon','common_name','scientifique_name','idg'});
ebd = unique(ebd,"sorted");

%% Match taxonomy

if false
    sp_ebird = groupcounts(ebd,["common_name","scientifique_name"]);
    
    % find ADU
    sp_2019 = readtable("data/2019_checklist_of_the_birds_of_kenya_short.csv", 'TextType', 'string');
    
    % Perform ADU adjustment to the 2019. Not exact match (in the futur) but
    % currently ok
    sp_2019.ADU(sp_2019.sort==45)=1254; % Elgon's Francolin -> Moorland Francolin
    sp_2019.ADU(sp_2019.sort==1004)=3198; % Common Scrub White-eye -> Abyssinian White-eye
    
    
    % Manual matching
    m={'Buteo buteo'         , 'Buteo buteo vulpinus/menetriesi';
    'Cisticola aberrans'     , 'Cisticola aberrans [emini Group]';
    'Cisticola lais'         , 'Cisticola lais distinctus';
    'Crithagra reichardi'    , 'Crithagra reichardi striatipectus';
    'Hyliota australis'      , 'Hyliota australis slatini';
    'Mirafra cheniana'       , 'Mirafra albicauda';
    'Oenanthe hispanica'     , 'Oenanthe hispanica melanoleuca';
    'Phalacrocorax carbo'    , 'Phalacrocorax carbo lucidus';
    'Pogoniulus bilineatus'  , 'Pogoniulus bilineatus leucolaimus';
    'Scleroptila psilolaema' , 'Scleroptila psilolaema elgonensis'};
    
    snm = sp_ebird.scientifique_name;
    for i_sp=1:size(m,1)
        tmp = strcmp(m{i_sp,1},sp_ebird.scientifique_name);
        snm(tmp) = repmat({m{i_sp,2}},sum(tmp),1);
    end
    
    % find match with table
    [~,tmp] = ismember(snm, sp_2019.Clements__scientific_name);
    
    % Add ADU to sp list
    sp_ebird.ADU(:) = nan;
    sp_ebird.ADU(tmp>0) = sp_2019.ADU(tmp(tmp>0));

    %writetable(sp_ebird,'data/eBird/sp_ebird.xlsx')
end


%% Read matching file
sp_ebird = readtable('data/eBird/sp_ebird.xlsx');
[~,tmp] = ismember(ebd.scientifique_name, sp_ebird.scientifique_name);
ebd.SEQ = sp_ebird.SEQ(tmp);

ebd(isnan(ebd.SEQ),:)=[];

%% Spatial grid

[~,id_sp] = ismember(ebd.SEQ, sp_old.SEQ);

[~,id_lat]=min((g.lat-ebd.lat).^2,[],2);
[~,id_lon]=min((g.lon-ebd.lon).^2,[],2);
map_ebird = false(size(map_old));
id = sub2ind(size(map_ebird),id_lat(id_sp>0),id_lon(id_sp>0),id_sp(id_sp>0));
map_ebird(id)=true;

%% Save
save('data/ebirdatlas.mat',"map_ebird")
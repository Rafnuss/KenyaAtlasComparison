
load('data/grid')
load('data/oldatlas')

%% Read and filter ebird data

ebd0 = readtable("data/eBird/ebd_KE_relOct-2023/ebd_KE_relOct-2023.txt",'TextType','string');
ebd = ebd0;
ebd = ebd(year(ebd.OBSERVATIONDATE)>=2009 & year(ebd.OBSERVATIONDATE)<=2023,:);
% NOT USED anymore ebd = ebd((ebd.CATEGORY=="species" | ebd.CATEGORY=="issf" | ebd.COMMONNAME=="Fischer's x Yellow-collared Lovebird (hybrid)") & year(ebd.OBSERVATIONDATE)>2000,:);

% Filter checklist too long or over too large distance (grid uncertain
numel(unique(ebd.SAMPLINGEVENTIDENTIFIER(ebd.EFFORTDISTANCEKM>30)))
ebd(ebd.EFFORTDISTANCEKM>30,:)=[]; % ebd.DURATIONMINUTES>24*60 & 

[~,id_lat]=min((g.lat-ebd.LATITUDE).^2,[],2);
[~,id_lon]=min((g.lon-ebd.LONGITUDE).^2,[],2);
ebd.idg = sub2ind(size(g.LAT),id_lat,id_lon);


%% Compute coverage map

ebd_checklist = groupsummary(ebd,{'SAMPLINGEVENTIDENTIFIER','ALLSPECIESREPORTED','DURATIONMINUTES','EFFORTDISTANCEKM','PROTOCOLTYPE','NUMBEROBSERVERS','idg'});
ebd_grid = groupsummary(ebd_checklist,"idg","sum","DURATIONMINUTES");

coverage_ebird=nan(size(g.LAT));
coverage_ebird(ebd_grid.idg) = ebd_grid.sum_DURATIONMINUTES/60;
coverage_ebird(isnan(coverage_ebird))=0;


%% Combine species/grid
ebd = table(ebd.LATITUDE,ebd.LONGITUDE,ebd.COMMONNAME,ebd.SCIENTIFICNAME, ebd.CATEGORY, ebd.idg,'VariableNames',{'lat','lon','common_name','scientific_name','category','idg'});
ebd = unique(ebd,"sorted");

%% explore
if false
    groupsummary(ebd,"PROTOCOLTYPE")
    
    numel(unique(ebd.COMMONNAME))
    
    tmp = groupsummary(ebd,"OBSERVERID");
    tmp = sortrows(tmp,"GroupCount","descend");
    height(tmp)
    sum(tmp.GroupCount>100)

    ebd_checklist = groupsummary(ebd,{'ALLSPECIESREPORTED','DURATIONMINUTES','EFFORTDISTANCEKM','LATITUDE','LONGITUDE','PROTOCOLTYPE','NUMBEROBSERVERS','idg'});
    
    sum(ebd_checklist.ALLSPECIESREPORTED)
    
    
    tmp2 = groupsummary(ebd_checklist(ebd_checklist.ALLSPECIESREPORTED==1,:),"idg","sum",["DURATIONMINUTES",'EFFORTDISTANCEKM','NUMBEROBSERVERS']);
    tmp=nan(size(g.LAT));
    tmp(tmp2.idg) = tmp2.sum_DURATIONMINUTES;
    
    figure; imagesc(g.lon,g.lat,tmp)
end


%% Match taxonomy

sp_ebird = readtable('data/eBird/sp_ebird.xlsx','TextType','string');

if false
    sp_ebird0 = groupcounts(ebd0,["COMMONNAME","SCIENTIFICNAME", "CATEGORY"]);
    

    [~,tmp] = ismember(sp_ebird0.SCIENTIFICNAME, sp_ebird.scientific_name);

    % Not matched
    sp_ebird0(~tmp,:)

    sp_ebird0(contains(sp_ebird0.COMMONNAME,"Shearwater"),:)

    snm = sp_ebird.scientific_name;
    for i_sp=1:size(m,1)
        tmp = strcmp(m{i_sp,1},sp_ebird.scientific_name);
        snm(tmp) = repmat({m{i_sp,2}},sum(tmp),1);
    end
    
    % find match with table
    [~,tmp] = ismember(snm, sp_2019.Clements__scientific_name);
    
    % Add ADU to sp list
    sp_ebird.ADU(:) = nan;
    sp_ebird.ADU(tmp>0) = sp_2019.ADU(tmp(tmp>0));

    %writetable(sp_ebird,'data/eBird/sp_ebird.xlsx')


    % all of sp_ebird are matched in sp_base
    A = ismember(sp_base.SEQ, sp_ebird.SEQ);
    sp_base(~A,:)

    A = ismember(sp_ebird.SEQ,sp_base.SEQ);
    sp_ebird(~A,:)
end

[Lia, Locb] = ismember(ebd.scientific_name, sp_ebird.scientific_name);
% All species should be match
tmp2 = unique(ebd((ebd.category=="species" | ebd.category=="issf") & ~Lia, ["common_name", "scientific_name"]));
assert(height(tmp2)==0)

% We should still have some slash left
unique(ebd(ebd.category=="slash" & ~Locb, ["common_name", "scientific_name"]))

% We only keep species for which we have a match
ebd = ebd(Lia,:);

% Add SEQ number
ebd.SEQ = sp_ebird.SEQ(Locb(Lia));

unique(ebd.common_name(ebd.SEQ==0))


%% Spatial grid
[~,id_sp] = ismember(ebd.SEQ, sp_base.SEQ);

[~,id_lat]=min((g.lat-ebd.lat).^2,[],2);
[~,id_lon]=min((g.lon-ebd.lon).^2,[],2);
map_ebird = false(size(map_old));
id = sub2ind(size(map_ebird),id_lat(id_sp>0),id_lon(id_sp>0),id_sp(id_sp>0));
map_ebird(id)=true;

%% Save
save('data/ebirdatlas.mat',"map_ebird", "coverage_ebird")
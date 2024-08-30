addpath(genpath("functions/"))
load('data/grid')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
load('data/grid_corr')

% Compute the difference
map_diff = (map_kbm|map_ebird)*2 - map_old; % convert to 2, 1, 0 and -1
map_diff(map_diff==0)=nan;
map_diff(map_diff==1)=0;
map_diff(map_diff==2)=1;

map_new_diff = map_kbm*2 - map_ebird; % convert to 2, 1, 0 and -1
map_new_diff(map_new_diff==0)=nan;
map_new_diff(map_new_diff==1)=0;
map_new_diff(map_new_diff==2)=1;

%% export map image
if true

szn = @(x, xs) sign(x).*min(sqrt(abs(x))/sqrt(abs(xs)),1);
szs = @(x, szmin, szmax) (szn(x,.7)+1)/2*(szmax-szmin)+szmin;
sz = @(x) szs(x, 50, 1000);

i_sp = find(strcmp(sp_base.common_name,"Great Crested Grebe"));
% i_sp = find(strcmp(sp_base.common_name,"House Sparrow"));

figure('position',[0 0 600 630]); 
tiledlayout('flow','Padding','none'); nexttile; hold on; box on; grid on

tmp = map_diff(:,:,i_sp); 
tmp2 = tmp; tmp2(tmp2==0)=1; tmp2(isnan(tmp2))=-1;

s0=scatter(g.LON(~g.mask),g.LAT(~g.mask),sz(tmp2(~g.mask).*g.corr(~g.mask)),'k','filled','markerfacealpha',0.3);

id =~isnan(tmp) & ~g.mask;
s1= scatter(g.LON(id), g.LAT(id),sz(tmp2(id).*g.corr(id).*~g.mask(id)),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);

id =~isnan(tmp) & g.mask;
s2= scatter(g.LON(id), g.LAT(id), sz(-3),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    

xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
axis equal tight; set(gca,"YDir","normal");
% title(sp_base.common_name(i_sp))
plot_google_map('MapType','terrain', 'ShowLabels',0);%[roadmap,,  satellite, terrain, hybrid);
colormap(brewermap(3,'RdYlGn'))
xticklabels([]); yticklabels([])
clim([-1 1])
axis([32.9394   42.4606   -4.9500    5.0500])

ann = annotation('textbox',Position=[0 0 .55 .08], String=sp_base.common_name(i_sp), Fontsize=16, FontWeight="bold", BackgroundColor='w');
ann2 = annotation('textbox',Position=[0 0 .55 .04], String=sp_base.scientific_name(i_sp), Fontsize=12, FontAngle='italic', EdgeColor="none");

% sp_base.IUCN(i_sp)
% iucn_match = [["Critically Endangered" "Data Deficient" "Endangered"  "Least Concern"  "Near Threatened"  "Vulnerable"]; ["CR" "DD" "EN" "LC" "NT" "VU"] ];
% tmp = iucn_match(2,iucn_match(1,:)==sp_base.IUCN(i_sp));
% FG = imread("export/iucn_"+tmp+".png");

% [~,id] = max(strlength(sp_base.scientific_name)) 25

for i_sp = 300:height(sp_base)
     
    delete(s0); delete(s1); delete(s2);

    tmp = map_diff(:,:,i_sp); 
    tmp2 = tmp; tmp2(tmp2==0)=1; tmp2(isnan(tmp2))=-1;

    s0=scatter(g.LON(~g.mask),g.LAT(~g.mask),sz(tmp2(~g.mask).*g.corr(~g.mask)),'k','filled','markerfacealpha',0.3);
    
    id =~isnan(tmp) & ~g.mask;
    s1= scatter(g.LON(id),g.LAT(id),sz(tmp2(id).*g.corr(id)),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    
    id =~isnan(tmp) & g.mask;
    s2= scatter(g.LON(id),g.LAT(id),sz(-3),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    
    ann.String=sp_base.common_name(i_sp);
    ann2.String=sp_base.scientific_name(i_sp);
    %title(sp_base.SEQ(i_sp) + " | "+sp_base.common_name(i_sp))
    exportgraphics(gca, "export/species/"+num2str(sp_base.SEQ(i_sp))+"_"+strrep(strrep(sp_base.common_name(i_sp),"/","_")," ","-")+".png")
    % keyboard
end
end




%% Method mfigure

% load data
ebd0 = readtable("data/eBird/ebd_KE_relOct-2023/ebd_KE_relOct-2023.txt",'TextType','string');
load('data/kbm/map_kbm')

% filter by years

% KBM
load('data/kbm/map_kbm')

%%
name = "White-headed vulture";


name = "African Rail";
name = "Crowned Crane";
name = "Senegal Plover";
name = "Common Pratincole";
name = validatestring(name,sp_base.common_name);

i_sp = find(strcmp(sp_base.common_name,name));

% filter species
ebd_sp = ebd0;
ebd_sp = ebd_sp(ebd_sp.COMMONNAME==sp_base.clements_common_name(i_sp),:);

ebd_y = ebd_sp;
ebd_y = ebd_y(year(ebd_y.OBSERVATIONDATE)>=2009 & year(ebd_y.OBSERVATIONDATE)<=2023,:);


i_adu = find(seq_order==sp_base.ADU(i_sp));

modes=["eBird", "eBird-upscaled", "combined", "KBM", "KBM-upscaled",  "old"];
figure('position',[0 0 1200 600]); tiledlayout('flow','Padding','none');

for i_m=1:numel(modes)
    mode = modes(i_m);
    ax(i_m) = nexttile;
hold on; box on; grid on
bck_col = [240 237 232]/255; % 
edg_col = [177 150 134]/255;
set(gca,'color',[137 192 229]/255);
% colormap(viridis)
list_country = ["Kenya", "United Republic of Tanzania", "Uganda", "Somalia", "Ethiopia", "Sudan"];
for i_c=1:numel(list_country)
    borders(list_country(i_c),'LineWidth',2,'facecolor',bck_col,'EdgeColor',edg_col)
end

if mode == "eBird"
    scatter(ebd_sp.LONGITUDE, ebd_sp.LATITUDE,15, [138 212 255]/255, 'o', 'filled')
    scatter(ebd_y.LONGITUDE, ebd_y.LATITUDE,15, [76 168 0]/255, 'o', 'filled')
elseif mode == "eBird-upscaled"
    imagesc(g.lon,g.lat,map_ebird(:,:,i_sp),'alphadata',.8*map_ebird(:,:,i_sp)); 
    colormap(ax(i_m), [76 168 0]/255)
elseif mode == "KBM"
    tmp = adhocp(:,:,i_adu);
    tmp(tmp>0) = 101;
    imagesc(gn.lon,gn.lat,tmp,'alphadata',adhocp(:,:,i_adu)>0)
    imagesc(gn.lon,gn.lat,fullpp(:,:,i_adu),'alphadata',fullp(:,:,i_adu)>0); 
    colormap(ax(i_m), [ ...
        repmat([254 246 114],25,1);
    repmat([254 226 114],25,1);
    repmat([254 205 114],50,1);
    repmat([254 167 114],100,1);
    repmat([254 151 89],100,1);
    repmat([255 128 51],200,1);
    repmat([255 68 0],250,1);
    repmat([208 0 98],250,1); ...
    repmat([255 161 195],1,1)]./[255 255 255])
    clim([0 101])
    colorbar
elseif mode == "KBM-upscaled"
    imagesc(g.lon,g.lat,map_kbm(:,:,i_sp),'alphadata',.8*map_kbm(:,:,i_sp)); 
    colormap(ax(i_m),[37 68 65]/255)
elseif mode == "combined" 
    imagesc(g.lon,g.lat,map_new_diff(:,:,i_sp),'alphadata',.8*(~isnan(map_new_diff(:,:,i_sp)))); 
    colormap(ax(i_m),[76 168 0; 50 110 30  ; 37 68 65]/255)
elseif mode == "old" 
    imagesc(g.lon,g.lat,map_old(:,:,i_sp),'alphadata',.8*map_old(:,:,i_sp)); 
    colormap(ax(i_m),[72, 50, 72]/255)
end

xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
xline([g.lon(1)-g.res/2 g.lon+g.res/2],Color=[.8 .8 .8], lineWidth=0.8)
yline([g.lat(1)-g.res/2 g.lat+g.res/2],Color=[.8 .8 .8], lineWidth=0.8)
set(gca,"YDir","normal");
axis equal tight;
axis([32.9394   42.4606   -4.9500    5.0500])
xticklabels([]); yticklabels([])

end

% title(sp_base.common_name(i_sp))

%% Difference

figure('position',[0 0 600 600]); hold on; box on; grid on
bck_col = [240 237 232]/255; % 
edg_col = [177 150 134]/255;
set(gca,'color',[137 192 229]/255);
colormap(viridis)
% borders("Kenya",'facecolor',bck_col,'EdgeColor',edg_col)
% borders("United Republic of Tanzania",'facecolor',bck_col,'EdgeColor',edg_col)
% borders("Uganda",'facecolor',bck_col,'EdgeColor',edg_col)
% borders("Somalia",'facecolor',bck_col,'EdgeColor',edg_col)
% borders("Ethiopia",'facecolor',bck_col,'EdgeColor',edg_col)
% borders("Sudan",'facecolor',bck_col,'EdgeColor',edg_col)
im = imagesc(g.lon,g.lat,map_diff(:,:,i_sp),'alphadata',.8*(~isnan(map_diff(:,:,i_sp)))); 
tmp = map_diff(:,:,i_sp); tmp(tmp==0)=1;tmp=tmp.*g.corr; id = ~isnan(tmp);
%scatter(g.LON(id),g.LAT(id),100,tmp(id)*8,'filled','MarkerEdgeColor','k')
%xline([g.lon(1)-g.res/2 g.lon+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
%yline([g.lat(1)-g.res/2 g.lat+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
axis equal tight; set(gca,"YDir","normal");
% axis([33.4500   41.9500   -5.5500    5.0500])
title(sp_base.common_name(i_sp))
plot_google_map('MapType','terrain', 'ShowLabels',0);%[roadmap,,  satellite, terrain, hybrid);
colormap(brewermap(3,'RdYlGn'))
xticklabels([]); yticklabels([])

%%

figure('position',[0 0 600 600]); hold on; box on; grid on
xline([g.lon(1)-g.res/2 g.lon+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
yline([g.lat(1)-g.res/2 g.lat+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
axis([32.9394   42.4606   -4.9500    5.0500])
plot_google_map('MapType','satellite', 'ShowLabels',0);%[roadmap,,  satellite, terrain, hybrid);
xticklabels([]); yticklabels([])

figure('position',[0 0 600 600]); hold on; box on; grid on
xline([g.lon(1)-g.res/2 g.lon+g.res/2],Color=[.8 .8 .8], LineWidth=2)
yline([g.lat(1)-g.res/2 g.lat+g.res/2],Color=[.8 .8 .8], LineWidth=2)
x=36.45; y=-1.45;
axis([x   x+3*.5   y   y+3*.5])
plot_google_map('MapType','satellite', 'ShowLabels',0);%[roadmap,,  satellite, terrain, hybrid);
xticklabels([]); yticklabels([])


%% 
%% Import data

load('data/grid')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")

%% Overview of the Old vs New atlas
map_sp_old = sum(map_old,3);
map_sp_new = sum(map_kbm|map_ebird,3);

figure('position',[0 0 1400 1200]); tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
colormap(viridis)
nexttile;hold on;
imagesc(g.lon,g.lat,map_sp_old,'alphadata',.8*(map_sp_old>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Species in old atlas')
plot_google_map('MapType','terrain', 'ShowLabels',0); box on; colorbar; caxis([0 630])

nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_kbm,3),'alphadata',.8*(sum(map_kbm,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Species in KBM')
plot_google_map('MapType','terrain', 'ShowLabels',0); box on; colorbar; caxis([0 630])

nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_ebird,3),'alphadata',.8*(sum(map_ebird,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Species in eBird')
plot_google_map('MapType','terrain', 'ShowLabels',0); box on; colorbar; caxis([0 630])

nexttile;hold on; title('Difference NEW (KBM+eBird) - OLD')
imagesc(g.lon,g.lat,map_sp_new-map_sp_old,'alphadata',.8*((map_sp_new-map_sp_old)~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar;
colormap(gca,brewermap([],'RdYlGn')); caxis([-300 300]);box on; 

nexttile;hold on; title('Difference KBM - OLD')
tmp = sum(map_kbm>0,3)-sum(map_old>0,3);
imagesc(g.lon,g.lat,tmp,'alphadata',.8*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar;
colormap(gca,brewermap([],'RdYlGn')); caxis([-300 300]); box on; 

nexttile;hold on; title('Difference eBird - OLD')
tmp = sum(map_ebird>0,3)-sum(map_old>0,3);
imagesc(g.lon,g.lat,tmp,'alphadata',.8*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar;
colormap(gca,brewermap([],'RdYlGn')); caxis([-300 300]); box on; 



%% Check map for a single species
i_sp = find(strcmp(sp_old.CommonName,"Black Kite"));

figure('position',[0 0 1400 600]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_old(:,:,i_sp),'alphadata',.8*map_old(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('Atlas')
plot_google_map('MapType','terrain', 'ShowLabels',0); 

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_kbm(:,:,i_sp),'alphadata',.8*map_kbm(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('KBM')
plot_google_map('MapType','terrain', 'ShowLabels',0);

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_ebird(:,:,i_sp),'alphadata',.8*map_ebird(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('eBird')
plot_google_map('MapType','terrain', 'ShowLabels',0);


%% Compute the difference
map_diff = (map_kbm|map_ebird)*2 - map_old; % convert to 2, 1, 0 and -1
map_diff(map_diff==0)=nan;
map_diff(map_diff==1)=0;
map_diff(map_diff==2)=1;
% nan: neither new nor old
% -1: only old
% 0: new and old
% 1: only new

%% Coverage
%opts = weboptions("Timeout",60);
% websave("data/birdmap/coverage.geojson","https://api.birdmap.africa/sabap2/v2/coverage/project/kenya/species?format=geojson",opts)
%cov_new = readtable("data/birdmap/coverage.csv", 'TextType', 'string');

% mask_table = table(g.SqN(:)+g.SqL(:), g.LON(:), g.LAT(:), reshape(sum(map_kbm|map_ebird,3),[],1), reshape(sum(map_old,3),[],1),...
%     VariableNames=["Sq","lon","lat","nb_kbmebird","nb_old"]);
% mask_table = mask_table(mask_table.Sq~="0",:);
% writetable(mask_table,"mask_table.xlsx")

% Histogram of species number
figure; hold on;
histogram(map_sp_old(map_sp_old>0),binWidth=10)
histogram(map_sp_new(map_sp_new>0),binWidth=10)

% Coverage & Effort
figure('position',[0 0 1400 600]); tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
colormap(viridis)
nexttile;hold on;
tmp =reshape(grp2idx(categorical(coverage_old(:))),size(coverage_old));
imagesc(g.lon,g.lat,tmp,'alphadata',0.8*(tmp>1)); 
axis equal tight; set(gca,"YDir","normal"); title('OLD: modeled coverage %')
plot_google_map('MapType','terrain', 'ShowLabels',0); c=colorbar;c.Ticks=0:6; c.TickLabels=unique(coverage_old)+"%";
xticks(''); yticks(''); box on
nexttile;hold on;
imagesc(g.lon,g.lat,log(coverage_kbm+coverage_ebird),'alphadata',.8*((coverage_ebird(:,:,1)+coverage_kbm(:,:,1))>0)); 
axis equal tight; set(gca,"YDir","normal"); title('NEW: sum of durations')
plot_google_map('MapType','terrain', 'ShowLabels',0); c=colorbar; 
c.Ticks=log([0 12 24 24*7 24*30.5 24*30.5*6]);
c.TickLabels=["0" "12hr" "1day" "1week" "1month" "6months"]; 
xticks(''); yticks(''); box on


% Fitted coverage
figure('position',[0 0 900 450]);tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile; hold on; title('Old Atlas')
x=categorical(coverage_old(:));
y=map_sp_old(:);
id = ~isundefined(x);
x=x(id);
y=y(id);
plot(x, y,'.k',MarkerSize=10)
p=fit(x, y, 'poly2');
plot(0:.1:7, p(0:.1:7),'-r')
xlabel("Coverage")
box on; grid on; ylim([0 700]); ylabel("Number of species")
% scatter(categorical(coverage_old(:)), reshape(sum(map_kbm|map_ebird,3),1,[]),'filled')

nexttile; hold on; title('New Atlas')
%plot(coverage_kbm(:), reshape(sum(map_kbm,3),[],1),'.b')
%plot(coverage_ebird(:), reshape(sum(map_ebird,3),[],1),'.r')
x=coverage_kbm(:)+coverage_ebird(:);
y=map_sp_new(:);
id = x~=0 & ~isnan(y);
x=x(id);
y=y(id);
plot(x, y,'.k',MarkerSize=10)
p=fit(log(x), y, 'poly2');
plot(logspace(0,5,100), p(log(logspace(0,5,100))),'-r')
set(gca,"XScale","log"); xlabel("Total duration");
set(gca,"Xtick", ([0 12 24 24*7 24*30.5 24*30.5*6]));
set(gca,"XtickLabels",["0" "12hr" "1day" "1week" "1month" "6months"]);
box on; grid on; ylim([0 700]); yticklabels('')


% Effort new vs old 
figure('position',[0 0 900 450]);
x=coverage_kbm(:)+coverage_ebird(:);
y=categorical(coverage_old(:));
% z = (reshape(sum(map_kbm|map_ebird,3),[],1) - reshape(sum(map_old,3),[],1)) ./ 2./(reshape(sum(map_kbm|map_ebird,3),[],1) + reshape(sum(map_old,3),[],1));
z = (reshape(sum(map_kbm|map_ebird,3),[],1) - reshape(sum(map_old,3),[],1));
id = ~isundefined(y) & y~="0";
x=x(id);
y=y(id);
z = z(id);
scatter(x, y, 300, z, 'filled',MarkerFaceAlpha=.8)
set(gca,"XScale","log"); xlabel("Total duration (new atlas)"); ylabel("Coverage (old atlas)")
set(gca,"Xtick", ([0 12 24 24*7 24*30.5 24*30.5*6]));
set(gca,"XtickLabels",["0" "12hr" "1day" "1week" "1month" "6months"]);
box on; grid on; colorbar; colormap(brewermap([],'RdYlGn'))
% clim([-.5 .5]); 
clim([-200 200])


% Compute correction 
% map_new_correction quantifies how many species in a square can be
% explained by the additional effort (compared to a standard effort for this category). 
figure('position',[0 0 900 400]); hold on;
cov_u = unique(coverage_old);
col = colormap( brewermap(7,'GnBu'));
map_new_correction = nan(size(map_kbm,1), size(map_kbm,2));
for i = 3:numel(cov_u)
    id = find(coverage_old(:)==cov_u(i));
    x = coverage_kbm(id)+coverage_ebird(id);
    id = id(x>0);
    x = x(x>0);
    y = sum(map_kbm|map_ebird,3) - sum(map_old,3);
    y=y(id);
    plot(x,y,'.', color=col(i,:), markersize=20)
    p = fit(log(x), y, 'poly1');
    x_interp = logspace(-0.6,4.5,100);
    plot(x_interp, p(log(x_interp)),'-', color=col(i,:), linewidth=2)
    p11 = predint(p,log(x_interp),0.95,'functional','on');
    fill([x_interp' ; flipud(x_interp')],[p11(:,1) ; flipud(p11(:,2))],col(i,:), EdgeColor = 'none',FaceAlpha=.2);  

    map_new_correction(id) = p(log(x));
end
set(gca,"XScale","log"); xlabel("Total duration"); ylabel("Number of species")
set(gca,"Xtick", ([0 12 24 24*7 24*30.5 24*30.5*6]));
set(gca,"XtickLabels",["0" "12hr" "1day" "1week" "1month" "6months"]);
yline(0,'--k')
box on; grid on; axis tight; ylim([-200 200])


% Corrected number of species
map_new_corrected = sum(map_kbm|map_ebird,3) - map_new_correction;

figure('position',[0 0 900 600]); tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
colormap( brewermap([],'RdYlGn')); 
nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_kbm|map_ebird,3)-sum(map_old,3),'alphadata',.8*(sum(map_kbm|map_ebird,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Old - uncorrected new')
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar; caxis([-200 200])

nexttile;hold on;
imagesc(g.lon,g.lat,map_new_corrected-sum(map_old,3),'alphadata',.8*(map_new_corrected>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Old - corrected new')
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar; caxis([-200 200])


%% Mask

mask = ~((coverage_kbm+coverage_ebird)>24) | coverage_old=="" | coverage_old=="0" | coverage_old=="11-30";
%mask_2 = ~((coverage_kbm+coverage_ebird)>24*7) | coverage_old=="" | coverage_old=="0" | coverage_old=="11-30" | coverage_old=="31-50";
%mask=zeros(size(mask_1));
%mask(~mask_1)=1;
%mask(~mask_2)=2;
% mask = sum(map_kbm|map_ebird,3)<100 | sum(map_old,3)<100;
% mask = abs(sum((map_kbm|map_ebird)>0,3)-sum(map_old>0,3))>100;

% figure;hold on; 
% imagesc(g.lon,g.lat,ones(size(mask)),'alphadata',.8*(1-mask/2)); 
% axis equal tight; set(gca,"YDir","normal")
% plot_google_map('MapType','terrain', 'ShowLabels',0); 
% set(gca,'XColor', 'none','YColor','none')
% set(gca, 'color', 'none');
% colormap( 'gray'); 


% Correction
% the correction C can be view as the probability that a given species is
% observed/not observed due to the change in effort. 
% The negative sign is to provide 

% C = - 10 * map_new_correction ./ (map_sp_new+map_sp_old)/2;
Clost = - map_new_correction ./ map_sp_new;
Cgain = - map_new_correction ./ map_sp_old;
C = - map_new_correction ./ (map_sp_new+map_sp_old)/2;


% Create size of marker for figure
sz = @(x) min(max(200*(1+4*x),15),600); % -1<x<1

figure('position',[0 0 600 600]);
colormap( brewermap([],'RdYlGn')); 
nexttile;hold on;
scatter(g.LON(~isnan(C)),g.LAT(~isnan(C)),sz(-C(~isnan(C))),'k','filled','MarkerEdgeColor','k')
% imagesc(g.lon,g.lat,nan(size(C)),'alphadata',zeros(size(C))); 
imagesc(g.lon,g.lat,ones(size(mask)),'alphadata',.6*mask); 
axis equal tight; set(gca,"YDir","normal");
plot_google_map('MapType','terrain', 'ShowLabels',0); colorbar; clim([-.1 .1])
box on; xticks(''); yticks('')
colormap(flipud(gray))


%%
i_sp = find(strcmp(sp_old.CommonName,"Great Crested Grebe"));
% i_sp = find(strcmp(sp_old.CommonName,"House Sparrow"));

figure('position',[0 0 600 600]); hold on; box on; grid on

tmp = map_diff(:,:,i_sp); 
tmp2 = tmp; tmp2(tmp2==0)=1; tmp2(isnan(tmp2))=-1;

s0=scatter(g.LON(~mask),g.LAT(~mask),sz(tmp2(~mask).*C(~mask)),'k','filled','markerfacealpha',0.3);

id =~isnan(tmp) & ~mask;
s1= scatter(g.LON(id),g.LAT(id),sz(tmp2(id).*C(id).*~mask(id)),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);

id =~isnan(tmp) & mask;
s2= scatter(g.LON(id),g.LAT(id),sz(-3),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    

xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
axis equal tight; set(gca,"YDir","normal");
title(sp_old.CommonName(i_sp))
plot_google_map('MapType','terrain', 'ShowLabels',0);%[roadmap,,  satellite, terrain, hybrid);
colormap(brewermap(3,'RdYlGn'))
xticklabels([]); yticklabels([])
clim([-1 1])
axis([32.9394   42.4606   -4.9500    5.0500])

for i_sp = 1:height(sp_old)
     
    delete(s0); delete(s1); delete(s2);

    tmp = map_diff(:,:,i_sp); 
    tmp2 = tmp; tmp2(tmp2==0)=1; tmp2(isnan(tmp2))=-1;

    s0=scatter(g.LON(~mask),g.LAT(~mask),sz(tmp2(~mask).*C(~mask)),'k','filled','markerfacealpha',0.3);
    
    id =~isnan(tmp) & ~mask;
    s1= scatter(g.LON(id),g.LAT(id),sz(tmp2(id).*C(id)),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    
    id =~isnan(tmp) & mask;
    s2= scatter(g.LON(id),g.LAT(id),sz(-3),tmp(id),'filled','MarkerEdgeColor',[.3 .3 .3]);
    

    title(sp_old.SEQ(i_sp) + " | "+sp_old.CommonName(i_sp))
    exportgraphics(gca, "export/species/"+num2str(sp_old.SEQ(i_sp))+"_"+strrep(sp_old.CommonName(i_sp)," ","-")+".png")
end

%%
% 
sp = sp_old;

tmp_old = map_old;
tmp_old(repmat(mask,1,1,size(tmp_old,3)))=0;
sp.old = squeeze(sum(tmp_old,[1 2]));
tmp_new = map_kbm|map_ebird;
tmp_new(repmat(mask,1,1,size(tmp_new,3)))=0;
sp.new = squeeze(sum(tmp_new,[1 2]));

sp.lost = squeeze(sum(tmp_old & ~tmp_new,[1 2]));
sp.kept = squeeze(sum(tmp_old & tmp_new,[1 2]));
sp.gain = squeeze(sum(~tmp_old & tmp_new,[1 2]));


% Convert C to a correction
% -1<C<1 -> 0<corr<1
% corr = max(min((C+1)/2,1),0);

tmp_lost = tmp_old & ~tmp_new;
tmp_gain = ~tmp_old & tmp_new;

sp.lost_score = squeeze(nansum(tmp_lost.*(1+C),[1 2]));
sp.gain_score = squeeze(nansum(tmp_gain.*(1+C),[1 2]));

% figure; hold on
% plot(sp.gain,sp.gain_score-sp.gain,'.g')
% plot(sp.lost,sp.lost_score-sp.lost,'.r')

writetable(sp,"export/sp_lost_kept_gain.csv")


g.corr = C;
g.mask = mask>=1;
g.coverage_new=coverage_kbm+coverage_ebird;
g.coverage_old=coverage_old;
save('data/grid_corr',"g")



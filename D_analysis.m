%% 
%% Import data

load('data/grid')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")

%% 
% 

figure('position',[0 0 1400 600]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
colormap(viridis)
nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_old>0,3),'alphadata',.8*(sum(map_old>0,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Total Species in Atlas')
plot_google_map; colorbar; caxis([0 630])

nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_kbm>0,3),'alphadata',.8*(sum(map_kbm>0,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Total Species in KBM')
plot_google_map; colorbar; caxis([0 630])

nexttile;hold on;
imagesc(g.lon,g.lat,sum(map_ebird>0,3),'alphadata',.8*(sum(map_ebird>0,3)>0)); 
axis equal tight; set(gca,"YDir","normal"); title('Total Species in eBird')
plot_google_map; colorbar; caxis([0 630])


figure('position',[0 0 1400 600]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
nexttile;hold on; title('Difference Atlas - KBM')
tmp = sum(map_old>0,3)-sum(map_kbm>0,3);
imagesc(g.lon,g.lat,tmp,'alphadata',.8*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map; colorbar;
colormap( brewermap([],'RdYlGn')); caxis([-300 300])

nexttile;hold on; title('Difference Atlas - eBird')
tmp = sum(map_old>0,3)-sum(map_ebird>0,3);
imagesc(g.lon,g.lat,tmp,'alphadata',.8*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map; colorbar;
colormap( brewermap([],'RdYlGn')); caxis([-300 300])

nexttile;hold on; title('Difference Atlas - (KBM+eBird)')
tmp = sum(map_old>0,3)-sum((map_kbm|map_ebird)>0,3);
imagesc(g.lon,g.lat,tmp,'alphadata',.6*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map; colorbar;
colormap( brewermap([],'RdYlGn')); caxis([-300 300])

%%
i_sp = find(strcmp(sp_old.CommonName,"Black Kite"));

figure('position',[0 0 1400 600]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
colormap(viridis)

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_old(:,:,i_sp),'alphadata',.8*map_old(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('Atlas')
plot_google_map; 

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_kbm(:,:,i_sp),'alphadata',.8*map_kbm(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('KBM')
plot_google_map;

nexttile;hold on; box on;
imagesc(g.lon,g.lat,map_ebird(:,:,i_sp),'alphadata',.8*map_ebird(:,:,i_sp)); 
axis equal tight; set(gca,"YDir","normal"); title('eBird')
plot_google_map;

%% Coverage
%opts = weboptions("Timeout",60);
% websave("data/birdmap/coverage.geojson","https://api.birdmap.africa/sabap2/v2/coverage/project/kenya/species?format=geojson",opts)
%cov_new = readtable("data/birdmap/coverage.csv", 'TextType', 'string');


%% Mask 
mask = sum(map_kbm|map_ebird,3)<100 | sum(map_old,3)<100;
% mask = abs(sum((map_kbm|map_ebird)>0,3)-sum(map_old>0,3))>100;

figure;hold on; title('Mask')
imagesc(g.lon,g.lat,ones(size(mask)),'alphadata',.8*(mask)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map; colormap( 'gray'); 

figure; hold on; title('Difference (new+eBird) - old')
tmp = sum((map_kbm|map_ebird)>0,3)-sum(map_old>0,3);
tmp(mask)=0;
imagesc(g.lon,g.lat,tmp,'alphadata',.8*(tmp~=0)); 
axis equal tight; set(gca,"YDir","normal")
plot_google_map; colorbar;
colormap( brewermap([],'RdYlGn')); caxis([-300 300])


%%
% 
diff = (map_kbm|map_ebird)*2 - map_old;
diff(repmat(mask,1,1,size(map_kbm,3)))=nan;
diff(diff==0)=nan;
diff(diff==1)=0;
diff(diff==2)=1;
% nan: neither new nor old
% -1: only old
% 0: new and old
% 1: only new


i_sp = find(strcmp(sp_old.CommonName,"Black Kite"));
figure('position',[0 0 600 600]); hold on; box on; grid on
colormap(viridis)
im = imagesc(g.lon,g.lat,diff(:,:,i_sp),'alphadata',.8*(~isnan(diff(:,:,i_sp)))); 
xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
axis equal tight; set(gca,"YDir","normal");
title(sp_old.CommonName(i_sp))
plot_google_map;
colormap(brewermap(3,'RdYlGn'))

for i_sp = 2:height(sp_old)
    im.CData = diff(:,:,i_sp);
    im.AlphaData = .8*(~isnan(diff(:,:,i_sp)));
    title(sp_old.SEQ(i_sp) + " | "+sp_old.CommonName(i_sp))
    exportgraphics(gca, "export/species/"+num2str(sp_old.SEQ(i_sp))+"_"+strrep(sp_old.CommonName(i_sp)," ","-")+".png")
end

%%
% 

sp_old.lost = squeeze(sum(diff==-1,[1 2]));
sp_old.kept = squeeze(sum(diff==0,[1 2]));
sp_old.gain = squeeze(sum(diff==1,[1 2]));
sp_old.diff = sp_old.gain-sp_old.lost;


writetable(sp_old,"export/lost_kept_gain_sp.csv")


figure; box on; grid on; grid on; grid on
histogram(sp_old.gain-sp_old.lost)
xlabel('Number of cell gained (+) or lost (-) since old atlas');
ylabel("Number of species")


sp_old_s=sortrows(sp_old,"diff");

figure; box on; grid on;
b = barh([sp_old.lost sp_old.gain].*[-1 1],'stacked');
c=brewermap(3,'RdYlGn');
b(1).FaceColor=  c(1,:);
b(2).FaceColor=  c(3,:);

nd=50;

for i_s=1:(nd*3):height(sp_old)
    figure('position',[0 0 3508 2480]/2); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    for u=1:3
        id_sub = (1:nd)+i_s-1+(u-1)*nd;
        id_sub = id_sub(id_sub<=height(sp_old));
        nexttile; box on; grid on;
        b = barh([-sp_old_s.kept(id_sub)/2 -sp_old_s.lost(id_sub) sp_old_s.kept(id_sub)/2 sp_old_s.gain(id_sub)],1,'stacked');
        c=brewermap(3,'RdYlGn');
        b(1).FaceColor=  c(2,:);
        b(2).FaceColor=  c(1,:);
        b(3).FaceColor=  c(2,:);
        b(4).FaceColor=  c(3,:);
        yticks(1:nd);
        yticklabels(sp_old_s.CommonName(id_sub));
        % xlabel('Number of grid lost (-) and gain(+)')
        set(gca, 'YDir','reverse')
        text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(sp_old_s.kept(id_sub)),'horiz','center'); 
        text(-sp_old_s.lost(id_sub)/2-sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.lost(id_sub)),'horiz','center'); 
        text(sp_old_s.gain(id_sub)/2+sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.gain(id_sub)),'horiz','center'); 
        grid on; axis tight; xlim([-90 90])
    end
    exportgraphics(gcf, "export/ranking_"+num2str(i_s)+".png")
    close gcf
end

%% Group species 

cat = {
'AM', 'Afrotropical migrant';...
'AMR', 'Afrotropical migrant and resident';...
'E', ' endemic species or race';...
%'EX', 'species which are thought to have become extinct in Kenya';...
%'HIST', 'no record for 50 years, i.e. no record since 1968 or earlier';...
%'IO', 'visitor from northwest Indian Ocean islands';...
%'MM', ' migrant from the Malagasy region';...
%'N', ' nomadic/wanderer';...
%'NRR', ' not recently recorded, i.e. since the period 1969 to 1999';...
%'OM', ' migrant from the Oriental region';...
'PM', ' migrant from the Palaearctic region';...
'PMR', ' migrant from the Palaearctic region and resident';...
%'RAR', ' Fewer than 5 records on East African Rarities Committee list at time of publication';...
%'RS', ' visitor from the Red Sea'...
%'SO', ' visitor from the Southern Ocean or Antarctica';...
%'VIO', ' vagrant from the northwest Indian Ocean';...
%'VM', ' vagrant from the Malagasy region';...
%'VN', ' vagrant from the Nearctic region';...
%'VO', ' vagrant from the Oriental region';...
%'VP', ' vagrant from the Palaearctic region';...
%'VSO', ' vagrant from the Southern Ocean or Antarctica';...
%'VSA', 'vagrant from southern Africa'...
};

for i_g=1:size(cat,1)
    id_sub = find(sp_old_s.(cat{i_g,1})=="TRUE");
    
    figure('position',[0 0 3508/3 2480]/2); tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp_old_s.kept(id_sub)/2 -sp_old_s.lost(id_sub) sp_old_s.kept(id_sub)/2 sp_old_s.gain(id_sub)],1,'stacked');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp_old_s.CommonName(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(sp_old_s.kept(id_sub)),'horiz','center'); 
    text(-sp_old_s.lost(id_sub)/2-sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.lost(id_sub)),'horiz','center'); 
    text(sp_old_s.gain(id_sub)/2+sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(cat{i_g,2})
    exportgraphics(gcf, "export/ranking_"+cat{i_g,1}+".png")
end

%% By status
cat = ["Critically Endangered","Endangered","Vulnerable","Near Threatened"];
for i_g=1:numel(cat)
    id_sub = find(sp_old_s.red_list==cat(i_g));
    
    figure('position',[0 0 3508/3 2480]/2); tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp_old_s.kept(id_sub)/2 -sp_old_s.lost(id_sub) sp_old_s.kept(id_sub)/2 sp_old_s.gain(id_sub)],1,'stacked');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp_old_s.CommonName(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(sp_old_s.kept(id_sub)),'horiz','center'); 
    text(-sp_old_s.lost(id_sub)/2-sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.lost(id_sub)),'horiz','center'); 
    text(sp_old_s.gain(id_sub)/2+sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(cat(i_g))
    exportgraphics(gcf, "export/ranking_"+cat(i_g)+".png")
end

%% By family

fam = groupcounts(sp_old_s,"family_english");
fam = fam(fam.GroupCount>10,:);

for i_f=1:height(fam)
    id_sub = find(sp_old_s.family_english==fam.family_english(i_f));
    
    figure('position',[0 0 3508/3 2480]/2); tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp_old_s.kept(id_sub)/2 -sp_old_s.lost(id_sub) sp_old_s.kept(id_sub)/2 sp_old_s.gain(id_sub)],1,'stacked');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp_old_s.CommonName(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(sp_old_s.kept(id_sub)),'horiz','center'); 
    text(-sp_old_s.lost(id_sub)/2-sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.lost(id_sub)),'horiz','center'); 
    text(sp_old_s.gain(id_sub)/2+sp_old_s.kept(id_sub)/2,1:numel(id_sub),num2str(sp_old_s.gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(fam.family_english(i_f))
    exportgraphics(gcf, "export/ranking_"+fam.family_english(i_f)+".png")
end

%%
G = groupsummary(sp_old_s,"family_english",@(x) median(x),["gain","lost","kept"]);
G.diff=G.fun1_gain-G.fun1_lost;
figure; hold on;
scatter(G.fun1_gain,G.fun1_lost,G.GroupCount*10,'filled')
text(G.fun1_gain,G.fun1_lost,G.family_english,"HorizontalAlignment","center")

G=sortrows(G,"diff");

figure('position',[0 0 3508/2 2480]/2); tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nd = 55;
for u=1:2
    id_sub = (1:nd)+(u-1)*nd;
    id_sub = id_sub(id_sub<=height(G));
    nexttile; box on; grid on;
    b = barh([-G.fun1_kept(id_sub)/2 -G.fun1_lost(id_sub) G.fun1_kept(id_sub)/2 G.fun1_gain(id_sub)],1,'stacked');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor =  c(2,:);
    b(2).FaceColor =  c(1,:);
    b(3).FaceColor =  c(2,:);
    b(4).FaceColor =  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(G.family_english(id_sub)+" ("+num2str(G.GroupCount(id_sub))+")");
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(G.fun1_kept(id_sub)),'horiz','center'); 
    text(-G.fun1_lost(id_sub)/2-G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_lost(id_sub)),'horiz','center'); 
    text(G.fun1_gain(id_sub)/2+G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-55 55])
end
exportgraphics(gcf, "export/ranking_family.png")
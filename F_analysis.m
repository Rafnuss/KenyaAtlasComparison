
%% Load data
addpath("functions/")
load('data/grid')
load('data/oldatlas')
load('data/kbmatlas.mat')
load("data/ebirdatlas.mat")
load('data/grid_corr')

c=brewermap(3,'RdYlGn');
c2=brewermap(11,'Paired');

%%

disp("Number of spcies/grid old atlas: "+ num2str(sum(map_old(:)>0)))
disp("Number of spcies/grid eBird: "+ num2str(sum(map_ebird(:)>0)))
disp("Number of spcies/grid KBM: "+ num2str(sum(map_kbm(:)>0)))
disp("Number of spcies/grid eBird+KBM: "+ num2str(sum(map_kbm(:)>0 | map_ebird(:)>0)))

disp("Number of spcies/grid old atlas: "+ num2str(sum(map_old>0 & ~g.mask,"all")))
disp("Number of spcies/grid eBird: "+ num2str(sum(map_ebird>0& ~g.mask,"all")))
disp("Number of spcies/grid KBM: "+ num2str(sum(map_kbm>0& ~g.mask,"all")))
disp("Number of spcies/grid eBird+KBM: "+ num2str(sum((map_kbm>0 | map_ebird>0)& ~g.mask,"all")))


%% Prepare data

sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");

% Make categorical value
sp.habitat = categorical(sp.habitat, ["Desert","Rock", "Grassland", "Shrubland", "Woodland", "Forest", "Human Modified", "Wetland", "Riverine", "Coastal", "Marine"]);
sp.trophic_niche = categorical(sp.trophic_niche, ["Frugivore", "Granivore", "Nectarivore", "Herbivore terrestrial", "Herbivore aquatic", "Invertivore", "Vertivore", "Aquatic predator", "Omnivore","Scavenger"]);
mig = categorical(["Sedentary", "Partial migratory", "Migratory"])';
sp.migration = mig(sp.migration);
hab_dens = categorical(["Dense habitats", "Semi-open habitats", "Open habitats"])';
sp.habitat_density = hab_dens(sp.habitat_density);

% Use the corrected value
% sp.gain = round(sp.gain_score);
% sp.lost = round(sp.lost_score);
% sp = removevars(sp,["lost_score","gain_score"]);

% Filteer for species with at least 10 squares total
% sp = sp((sp.lost+sp.gain+sp.kept)>10,:);

% Remove species which have not been recorded in both atlas
sp = sp(~(sp.kept + sp.gain + sp.lost)==0,:);

% Compute diff and sort
sp.diff = sp.new - sp.old;
% sp_s.diff_prop = sp_s.gain./(sp_s.new) -sp_s.lost./(sp_s.old);
sp.diff_2 = (sp.new - sp.old) ./ (sp.new + sp.old);
sp.diff_prop = (sp.new - sp.old) ./ sp.old;
sp.diff_prop_2 = (sp.new - sp.old) ./ (sp.kept + sp.gain + sp.lost);
sp = sortrows(sp,"diff");


%% All Species

nd=54;

sp=sortrows(sp,"diff");
% sp=sortrows(sp,"SEQ");

% height(sp_s)/3/nd

for i_s=1:(nd*3):height(sp)
    figure('position',[0 0 3508 2480]/2); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    for u=1:3
        id_sub = (1:nd)+i_s-1+(u-1)*nd;
        id_sub = id_sub(id_sub<=height(sp));
        nexttile; box on; grid on;
        b = barh([-sp.kept(id_sub)/2 -sp.lost(id_sub) sp.kept(id_sub)/2 sp.gain(id_sub)],1,'stacked',EdgeColor='none');
        c=brewermap(3,'RdYlGn');
        b(1).FaceColor = c(2,:);
        b(2).FaceColor = c(1,:);
        b(3).FaceColor = c(2,:);
        b(4).FaceColor = c(3,:);
        yticks(1:nd);
        yticklabels(sp.common_name(id_sub));
        % xlabel('Number of grid lost (-) and gainlost(+)')
        set(gca, 'YDir','reverse')
        text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp.kept(id_sub))),'horiz','center'); 
        text(-sp.lost(id_sub)/2-sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.lost(id_sub))),'horiz','center'); 
        text(sp.gain(id_sub)/2+sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.gain(id_sub))),'horiz','center'); 
        grid on; axis tight; xlim([-90 90])
    end
    %exportgraphics(gcf, "export/ranking/"+num2str(i_s)+".png")
    close gcf
end




%% Habitat, Trophic_Niche, Migration

c2=cell(3,1);
c2{1} = [
  [240, 230, 140],   % Desert
  [128, 128, 128],   % Rock
  [124, 252, 0],     % Grassland
  [85, 107, 47],     % Shrubland
  [139, 69, 19],     % Woodland
  [0, 100, 0],       % Forest
  [169, 169, 169],   % Human Modified
  [143, 151, 121],   % Wetland
  [65, 105, 225],    % Riverine
  [30, 144, 255],    % Coastal
  [0, 0, 128]        % Marine
]/255;

c2{2} =[
  [241, 196, 15],    % Frugivore
  [153, 102, 0],     % Granivore
  [255, 195, 0],     % Nectarivore
  [46, 204, 113],    % Herbivore terrestrial
  [0, 255, 255],     % Herbivore aquatic
  [52, 73, 94],      % Invertivore
  [155, 89, 182],    % Vertivore
  [0, 191, 255],     % Aquatic predator
  [255, 165, 0],     % Omnivore
  [128, 128, 128]    % Scavenger
]/255;

c2{3} =[
  [92, 53, 102],     % Sedentary
  [46, 134, 193],    % Partial migratory
  [244, 81, 30]      % Migratory
]/255;

c2{4} =[
  [139, 69, 19],     % dense
  [85, 107, 47],    % semi-open 
  [210, 180, 140]      % open
]/255;

c2{5} =[
  [135, 206, 250],     % Aerial
  [36, 157, 159],    % Aquatic
  [255, 215, 0],      % Generalist
  [210, 180, 140],    % Insessorial
  [139, 69, 19]      % Terrestrial
]/255;

c2{6} =[
  [153, 17, 17],     % Carnivore
  [46, 204, 113],    % Herbivore
  [255, 165, 0],     % Omnivore
  [128, 128, 128]    % Scavenger
]/255;


cat_list = ["habitat", "trophic_niche", "migration", "habitat_density", "primary_lifestyle", "trophic_level"];
for u = 1:numel(cat_list)
    figure('position',[0 0 600 (length(unique(sp.(cat_list(u))))+1)*40]);
    tt = char(strrep(cat_list(u), "_", " "));
    tt(1)=upper(tt(1));
    title(tt)
    violinplot(sp.diff,sp.(cat_list(u)), 'ViolinColor',c2{u}); view([90 -90])
    grid on; axis tight; ylim([-30 30]);
    set(gca,'xdir','reverse')
    set(gcf,'renderer','Painters')
    ylabel("Change in number of squares")
    exportgraphics(gcf, "export/traits/cat_"+cat_list(u)+".eps")
    exportgraphics(gcf, "export/traits/cat_"+cat_list(u)+".png")
end

figure; 
scatter(log10(sp.mass),sp.diff,'k','filled',"MarkerFaceAlpha",.2)
l=lsline;l.Color="r";

% figure; 
% scatter((sp.Range_Size),sp.diff,'k','filled',"MarkerFaceAlpha",.2)
% l=lsline;l.Color="r";

%%  Migrant
sum(sp.palearctic==1)
1 - sum(sp.new(sp.palearctic==1))/sum(sp.old(sp.palearctic==1))

mean(sp.diff(sp.palearctic==1))
mean(sp.diff(sp.palearctic~=1))


[h,p,ks2stat] = kstest2(sp.diff(sp.palearctic==1),sp.diff(sp.palearctic~=1),'tail','larger');

% Using avonet classification
sp(sp.palearctic==1 & sp.migration~="Migratory",:)
sp(sp.palearctic~=1 & sp.migration=="Migratory",:)
mean(sp.diff(sp.migration=="Migratory"))
mean(sp.diff(sp.migration~="Migratory"))

mean(sp.diff_prop(sp.migration=="Migratory"))


tmp = countlabels(sp.checklist_family(sp.palearctic==1));
fam = tmp.Label(tmp.Count>=2);

tmp2 = sp(ismember(sp.checklist_family,string(fam)),:);
tmp3 = groupsummary(tmp2, ["checklist_family", "palearctic"], ["mean", "std"], ["lost" "kept" "gain" "diff"]);

tmp32 = groupsummary(tmp3,"checklist_family","std","mean_diff");
tmp3 = join(tmp3,tmp32(:,["checklist_family", "std_mean_diff"]));
tmp3 = sortrows(tmp3,"std_mean_diff");

[~, tmp4] = ismember(tmp2.checklist_family, unique(tmp3.checklist_family));

figure; hold on;
b=barh(reshape(tmp3.mean_diff, 2, [])');
b(2).FaceColor = [252, 182, 68]/255;
b(1).FaceColor = [186, 120, 234]/255;
errorbar(tmp3.mean_diff, repelem(1:height(tmp3)/2,1,2)'+repmat([-.15 .15],1,height(tmp3)/2)',tmp3.std_diff,'horizontal',"LineStyle","none", "CapSize",3,"Color",[.7 .7 .7])
% scatter(tmp2.diff, tmp4, [], tmp2.palearctic, "filled")
% scatter(tmp3.mean_diff, repelem(1:10,1,2)', 100, tmp3.palearctic, "filled")
yticks(1:height(tmp3)/2)
yticklabels(unique(tmp3.checklist_family,"stable"))
 grid on; box on; xlim([-40 40]); set(gca, 'YDir','reverse');

exportgraphics(gcf, "figures/migrant_per_family.eps") 

% 
% figure; yline(0)
% violinplot(tmp2.diff,tmp2.checklist_family + tmp2.palearctic); view([90 -90])
% set(gca, 'XDir','reverse'); grid on; box on; axis tight;ylim([-50 50]); 
%% Group species 

cat = ["AfrotropicalMigrant" "Endemic" "Palearctic" "waterBird"];

for i_g=1:numel(cat)
    id_sub = find(sp.(cat(i_g)));
    nd = ceil(numel(id_sub)/3);
    figure('position',[0 0 3508 2480]/2); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    for u=1:3
        id_sub2 = (1:nd)+1-1+(u-1)*nd;
        id_sub2 = id_sub2(id_sub2<=numel(id_sub));
        id_sub2 = id_sub(id_sub2);

        nexttile; box on; grid on;
        b = barh([-sp.kept(id_sub2)/2 -sp.lost(id_sub2) sp.kept(id_sub2)/2 sp.gain(id_sub2)],1,'stacked',EdgeColor='none');
        c=brewermap(3,'RdYlGn');
        b(1).FaceColor=  c(2,:);
        b(2).FaceColor=  c(1,:);
        b(3).FaceColor=  c(2,:);
        b(4).FaceColor=  c(3,:);
        yticks(1:numel(id_sub2));
        yticklabels(sp.common_name(id_sub2));
        % xlabel('Number of grid lost (-) and gain(+)')
        set(gca, 'YDir','reverse')
        text(zeros(numel(id_sub2),1),1:numel(id_sub2),num2str(round(sp.kept(id_sub2))),'horiz','center'); 
        text(-sp.lost(id_sub2)/2-sp.kept(id_sub2)/2,1:numel(id_sub2),num2str(round(sp.lost(id_sub2))),'horiz','center'); 
        text(sp.gain(id_sub2)/2+sp.kept(id_sub2)/2,1:numel(id_sub2),num2str(round(sp.gain(id_sub2))),'horiz','center'); 
        grid on; axis tight; xlim([-90 90])
       % title(cat(i_g))
    end
    %exportgraphics(gcf, "export/status/"+cat(i_g)+".png")
end






%% By status
cat = ["Critically Endangered","Endangered","Vulnerable","Near Threatened"];
for i_g=1:numel(cat)
    id_sub = find(sp.IUCN==cat(i_g));
    
    figure('position',[0 0 3508/3 100+2480*numel(id_sub)/100]/2);  tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp.kept(id_sub)/2 -sp.lost(id_sub) sp.kept(id_sub)/2 sp.gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp.common_name(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp.kept(id_sub))),'horiz','center'); 
    text(-sp.lost(id_sub)/2-sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.lost(id_sub))),'horiz','center'); 
    text(sp.gain(id_sub)/2+sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.gain(id_sub))),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(cat(i_g))
    exportgraphics(gcf, "export/status/"+cat(i_g)+".png")
end

%% By family

fam = groupcounts(sp,"Family");
fam = fam(fam.GroupCount>1,:);

for i_f=1:height(fam)
    id_sub = find(sp.Family==fam.Family(i_f));
    
    figure('position',[0 0 3508/3 100+2480*numel(id_sub)/100]/2); tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp.kept(id_sub)/2 -sp.lost(id_sub) sp.kept(id_sub)/2 sp.gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp.common_name(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp.kept(id_sub))),'horiz','center'); 
    text(-sp.lost(id_sub)/2-sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.lost(id_sub))),'horiz','center'); 
    text(sp.gain(id_sub)/2+sp.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp.gain(id_sub))),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(fam.Family(i_f))
    exportgraphics(gcf, "export/family/"+fam.Family(i_f)+".eps")
end


%%
G = groupsummary(sp,"checklist_family",@(x) median(x),["gain","lost","kept"]);
G.diff_score=G.fun1_gain-G.fun1_lost;

figure; hold on; grid on; box on; axis equal
plot([0 100],[0 100],'-k')
scatter(G.fun1_kept+G.fun1_gain,G.fun1_kept+G.fun1_lost,G.GroupCount*10,'filled')
text(G.fun1_kept+G.fun1_gain,G.fun1_kept+G.fun1_lost,G.checklist_family,"HorizontalAlignment","center")

%% 
G=sortrows(G,"diff_score");

figure('position',[0 0 3508/2 2480]/2); tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nd = 46;
for u=1:2
    id_sub = (1:nd)+(u-1)*nd;
    id_sub = id_sub(id_sub<=height(G));
    nexttile; box on; grid on;
    b = barh([-G.fun1_kept(id_sub)/2 -G.fun1_lost(id_sub) G.fun1_kept(id_sub)/2 G.fun1_gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor =  c(2,:);
    b(2).FaceColor =  c(1,:);
    b(3).FaceColor =  c(2,:);
    b(4).FaceColor =  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(G.checklist_family(id_sub)+" ("+num2str(G.GroupCount(id_sub))+")");
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(G.fun1_kept(id_sub)),'horiz','center'); 
    text(-G.fun1_lost(id_sub)/2-G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_lost(id_sub)),'horiz','center'); 
    text(G.fun1_gain(id_sub)/2+G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-55 55])
end
exportgraphics(gcf, "export/family/all.png")


%% Scavanger

sp(sp.trophic_niche=="Scavenger"| sp.common_name=="Bateleur",["common_name", "diff_prop"])
% missing bateleur

sp(sp.trophic_level=="Scavenger"| sp.common_name=="Bateleur",["common_name", "diff_prop"])
% missing egyptian vulture and bateleur

%% Introduced

mean(sp.diff(sp.habitat=="Human Modified"))
sp.common_name(sp.habitat=="Human Modified")

sp_introduced = ["Fischer's x Yellow-collared Lovebird", "House Sparrow", "Indian House Crow", "Feral Pigeon"];

id = find(ismember(sp.common_name, sp_introduced));

[sp.common_name(id) sp.old(id) sp.new(id)]




%% Waterbird

sp_water = sp(sp.waterBird==1,:);

figure('position',[0 0 900 900]); hold on; grid on; box on; axis equal square
plot([0 100],[0 100],'--k')
scatter(sp_water.old, sp_water.new,100, categorical(sp_water.Family),'filled','MarkerFaceAlpha',.8,'MarkerEdgeColor','none')
xlabel('Number of square in old atlas'); ylabel('Number of square in new atlas')


%% General trend NOT USED 

% wider range -> wider increase
x = sp.gain +sp.lost+ sp.kept;
y = sp.gain - sp.lost;

p = fit(x, y, 'poly1');
x_interp = 0:100;
p11 = predint(p,x_interp,0.95,'functional','on');


figure; hold on;
fill([x_interp' ; flipud(x_interp')],[p11(:,1) ; flipud(p11(:,2))],col(i,:), EdgeColor = 'none',FaceAlpha=.2);  

plot(x_interp, p(x_interp),'-', linewidth=2)
scatter( x ,  y, 'ok','filled',AlphaData=0.5 )
yline(mean(sp.gain - sp.lost))
yline(0)
ylim([-25 25])

[rho, p]=corrcoef(sp.gain +sp.lost+ sp.kept , sp.gain - sp.lost);


% Globale histogram
figure; box on; grid on; grid on; grid on
histogram(sp.gain-sp.lost)
xlabel('Number of cell gained (+) or lost (-) since old atlas');
ylabel("Number of species")
% exportgraphics(gcf, "figures/histogram.png")


%
figure;  box on; grid on;
b = barh([sp.lost sp.gain].*[-1 1],'stacked');
b(1).FaceColor = c(1,:);
b(2).FaceColor = c(3,:);



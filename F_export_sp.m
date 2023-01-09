c=brewermap(3,'RdYlGn');

sp_s = readtable("export/sp_lost_kept_gain.csv", TextType="string");


% Use the corrected value
%sp_s.gain = round(sp_s.gain_score);
%sp_s.lost = round(sp_s.lost_score);
sp_s = removevars(sp_s,["lost_score","gain_score"]);

% Filteer for species with at least 10 squares total
sp_s = sp_s((sp_s.lost+sp_s.gain+sp_s.kept)>10,:);

% Compute diff and sort
sp_s.diff = sp_s.new -sp_s.old;
% sp_s.diff_prop = sp_s.gain./(sp_s.new) -sp_s.lost./(sp_s.old);
% sp_s.diff = (sp_s.new -sp_s.old) ./ (sp_s.new + sp_s.old);
sp_s=sortrows(sp_s,"diff");


%%
figure; box on; grid on; grid on; grid on
histogram(sp.gain-sp.lost)
xlabel('Number of cell gained (+) or lost (-) since old atlas');
ylabel("Number of species")


%% 
figure;  box on; grid on;
b = barh([sp.lost sp.gain].*[-1 1],'stacked');
b(1).FaceColor = c(1,:);
b(2).FaceColor = c(3,:);


%%

nd=54;

% height(sp_s)/3/nd

for i_s=1:(nd*3):height(sp_s)
    figure('position',[0 0 3508 2480]/2); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    for u=1:3
        id_sub = (1:nd)+i_s-1+(u-1)*nd;
        id_sub = id_sub(id_sub<=height(sp_s));
        nexttile; box on; grid on;
        b = barh([-sp_s.kept(id_sub)/2 -sp_s.lost(id_sub) sp_s.kept(id_sub)/2 sp_s.gain(id_sub)],1,'stacked',EdgeColor='none');
        c=brewermap(3,'RdYlGn');
        b(1).FaceColor = c(2,:);
        b(2).FaceColor = c(1,:);
        b(3).FaceColor = c(2,:);
        b(4).FaceColor = c(3,:);
        yticks(1:nd);
        yticklabels(sp_s.CommonName(id_sub));
        % xlabel('Number of grid lost (-) and gainlost(+)')
        set(gca, 'YDir','reverse')
        text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp_s.kept(id_sub))),'horiz','center'); 
        text(-sp_s.lost(id_sub)/2-sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.lost(id_sub))),'horiz','center'); 
        text(sp_s.gain(id_sub)/2+sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.gain(id_sub))),'horiz','center'); 
        grid on; axis tight; xlim([-90 90])
    end
    exportgraphics(gcf, "export/ranking/"+num2str(i_s)+".png")
    close gcf
end

%% Group species 
cat = ["AfrotropicalMigrant" "Endemic" "Palearctic"];

for i_g=1:numel(cat)
    id_sub = find(sp_s.(cat(i_g)));
    nd = ceil(numel(id_sub)/3);
    figure('position',[0 0 3508 2480]/2); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    for u=1:3
        id_sub2 = (1:nd)+1-1+(u-1)*nd;
        id_sub2 = id_sub2(id_sub2<=numel(id_sub));
        id_sub2 = id_sub(id_sub2);

        nexttile; box on; grid on;
        b = barh([-sp_s.kept(id_sub2)/2 -sp_s.lost(id_sub2) sp_s.kept(id_sub2)/2 sp_s.gain(id_sub2)],1,'stacked',EdgeColor='none');
        c=brewermap(3,'RdYlGn');
        b(1).FaceColor=  c(2,:);
        b(2).FaceColor=  c(1,:);
        b(3).FaceColor=  c(2,:);
        b(4).FaceColor=  c(3,:);
        yticks(1:numel(id_sub2));
        yticklabels(sp_s.CommonName(id_sub2));
        % xlabel('Number of grid lost (-) and gain(+)')
        set(gca, 'YDir','reverse')
        text(zeros(numel(id_sub2),1),1:numel(id_sub2),num2str(round(sp_s.kept(id_sub2))),'horiz','center'); 
        text(-sp_s.lost(id_sub2)/2-sp_s.kept(id_sub2)/2,1:numel(id_sub2),num2str(round(sp_s.lost(id_sub2))),'horiz','center'); 
        text(sp_s.gain(id_sub2)/2+sp_s.kept(id_sub2)/2,1:numel(id_sub2),num2str(round(sp_s.gain(id_sub2))),'horiz','center'); 
        grid on; axis tight; xlim([-90 90])
       % title(cat(i_g))
    end
    exportgraphics(gcf, "export/status/"+cat(i_g)+".png")
end

%% By status
cat = ["Critically Endangered","Endangered","Vulnerable","Near Threatened"];
for i_g=1:numel(cat)
    id_sub = find(sp_s.IUCN==cat(i_g));
    
    figure('position',[0 0 3508/3 100+2480*numel(id_sub)/100]/2);  tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp_s.kept(id_sub)/2 -sp_s.lost(id_sub) sp_s.kept(id_sub)/2 sp_s.gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp_s.CommonName(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp_s.kept(id_sub))),'horiz','center'); 
    text(-sp_s.lost(id_sub)/2-sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.lost(id_sub))),'horiz','center'); 
    text(sp_s.gain(id_sub)/2+sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.gain(id_sub))),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(cat(i_g))
    exportgraphics(gcf, "export/status/"+cat(i_g)+".png")
end

%% By family

fam = groupcounts(sp_s,"Family");
fam = fam(fam.GroupCount>=5,:);

for i_f=1:height(fam)
    id_sub = find(sp_s.Family==fam.Family(i_f));
    
    figure('position',[0 0 3508/3 100+2480*numel(id_sub)/100]/2); tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
    nexttile; box on; grid on;
    b = barh([-sp_s.kept(id_sub)/2 -sp_s.lost(id_sub) sp_s.kept(id_sub)/2 sp_s.gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:numel(id_sub));
    yticklabels(sp_s.CommonName(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(sp_s.kept(id_sub))),'horiz','center'); 
    text(-sp_s.lost(id_sub)/2-sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.lost(id_sub))),'horiz','center'); 
    text(sp_s.gain(id_sub)/2+sp_s.kept(id_sub)/2,1:numel(id_sub),num2str(round(sp_s.gain(id_sub))),'horiz','center'); 
    grid on; axis tight; xlim([-90 90])
    title(fam.Family(i_f))
    exportgraphics(gcf, "export/family/"+fam.Family(i_f)+".png")
end



%% 
G=sortrows(G,"diff_score");

figure('position',[0 0 3508/2 2480]/2); tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nd = 55;
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
    yticklabels(G.Family(id_sub)+" ("+num2str(G.GroupCount(id_sub))+")");
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(G.fun1_kept(id_sub)),'horiz','center'); 
    text(-G.fun1_lost(id_sub)/2-G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_lost(id_sub)),'horiz','center'); 
    text(G.fun1_gain(id_sub)/2+G.fun1_kept(id_sub)/2,1:numel(id_sub),num2str(G.fun1_gain(id_sub)),'horiz','center'); 
    grid on; axis tight; xlim([-55 55])
end
exportgraphics(gcf, "export/family/all.png")


%%
G = groupsummary(sp_s,"Family",@(x) median(x),["gain","lost","kept"]);
G.diff_score=G.fun1_gain-G.fun1_lost;

figure; hold on; grid on; box on; axis equal
plot([0 100],[0 100],'-k')
scatter(G.fun1_kept+G.fun1_gain,G.fun1_kept+G.fun1_lost,G.GroupCount*10,'filled')
text(G.fun1_kept+G.fun1_gain,G.fun1_kept+G.fun1_lost,G.Family,"HorizontalAlignment","center")

%%
figure('position',[0 0 900 900]); hold on; grid on; box on; axis equal square
plot([0 100],[0 100],'--k')
scatter(sp_s.old,sp_s.new,'k','filled','MarkerFaceAlpha',.4,'MarkerEdgeColor','none')
xlabel('Number of square in old atlas'); ylabel('Number of square in new atlas')

fam = "Hawks, Vultures, Buzzards, Eagles and Allies"; %Bee-eaters, Bustards Ducks and Geese
scatter(sp_s.old(sp_s.Family==fam),sp_s.new(sp_s.Family==fam),300,'r','filled','MarkerFaceAlpha',.6,'MarkerEdgeColor','none')
% text(sp_s.old(sp_s.Family==fam),sp_s.new(sp_s.Family==fam)-1.5,sp_s.CommonName(sp_s.Family==fam),HorizontalAlignment="center", VerticalAlignment="top", FontSize=18)




figure('position',[0 0 900 900]); hold on; grid on; box on; axis equal square
fam = "Hawks, Vultures, Buzzards, Eagles and Allies"; %Bee-eaters, Bustards Ducks and Geese
tmp = sortrows(sp_s(sp_s.Family==fam,:),"diff");
b=barh(max(tmp.diff,0),1); b.FaceColor=c(3,:);
b=barh(min(tmp.diff,0),1); b.FaceColor=c(1,:);
xlim([-1 1])
yticks(1:height(tmp)); yticklabels(tmp.CommonName)
% text(sp_s.old(sp_s.Family==fam),sp_s.new(sp_s.Family==fam)-1.5,sp_s.CommonName(sp_s.Family==fam),HorizontalAlignment="center", VerticalAlignment="top", FontSize=18)






%% Compare raptor
D = table(repelem([99,91,89,92,87,129,88]',4,1),...
    repmat(["1980";"1980";"2000";"2000"],7,1),...    
    repmat(["non-migration";"migration"],7*2,1),...
    [3.04  6.84 4.1 7.78 0.4 7.23 0 0 24.14 116.76 21.27 46.1 2.17 7.21 2.30 1.26 4.96 14.16 4.44 5.23 6.19 8.3 5.93 4.43 1.83 1.87 1.36 0.7]',...
    VariableNames=["SEQ" "years" "season","v"]);

D_comp = unstack(D,"v","years");
D_comp.diff = (D_comp.x2000 - D_comp.x1980) ./ D_comp.x1980;
D_comp = unstack(removevars(D_comp,["x1980", "x2000"]),"diff","season");

[~,id] = ismember(D_comp.SEQ,sp_s.SEQ);

figure; hold on;
bar(categorical(sp_s.CommonName(id(id>0))), [D_comp.migration D_comp.non_migration])

sp_vult = sp_s(ismember(sp_s.SEQ,[D.SEQ; 90]),:);
(sp_vult.new-sp_vult.old) ./sp_vult.old
sp_vult.CommonName


a=sp_s(sp_s.CommonName=="House Sparrow",:);
(a.new-a.old) ./a.old


%% On a tree
Tree = phytreeread("data/addTaxaComplete2Sep2022.tre");
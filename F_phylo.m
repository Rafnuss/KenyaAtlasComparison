

sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");
Tree = phytreeread("data/phylo/kenyan_bird.tre");

D = pdist(Tree, Squareform=1);

leaf = str2double(string(get(Tree, 'LeafNames')));

[~,id] = ismember(leaf,sp.SEQ);
sp_o = sp(id,:);


% D2 = squareform(pdist(sp_o.gain_score-sp_o.lost_score));
d = (sp_o.gain_score-sp_o.lost_score);
d = (d-mean(d))./std(d);
C = (d .* d');

figure; tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
D2=D;D2(D2>50)=nan; 
nexttile; imagesc(D2); axis equal square;
nexttile; imagesc(C); axis equal square;


dh=5;
dh2=1;
h=1:dh2:150;
v=nan(size(h));
for i_h=1:numel(h)
    id = D>=(h(i_h)-dh/2) & D<(h(i_h)+dh/2) & D~=0;
    v(i_h) = mean(C(id));
end

figure; hold on; box on; grid on;
plot(D(:),C(:),'.k')
plot(h,v,'-or',LineWidth=3)
xlim([1 150])
xlim([1 50])
ylim([-.1 .3])

%%

sp.diff = sp.new - sp.old;
sps2 = sp(ismissing(sp.flag),:); 
[~,id_min] = mink(sps2.diff,10);
[~,id_max] = maxk(sps2.diff,10);



figure('position',[0 0 600 200]); tiledlayout('flow','Padding','none');
for i=1:2
    if i==1
        sps = sps2(id_min,:);
        x_lim = [-80 50];
    else
        sps = sps2(id_max,:);
        x_lim = [-40 100];
    end
    nexttile;
    b = barh([-sps.kept(:)/2 -sps.lost(:) sps.kept(:)/2 sps.gain(:)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor = c(2,:);
    b(2).FaceColor = c(1,:);
    b(3).FaceColor = c(2,:);
    b(4).FaceColor = c(3,:);
    yticks(1:height(sps));
    yticklabels(sps.common_name(:));
    % xlabel('Number of grid lost (-) and gainlost(+)')
    set(gca, 'YDir','reverse')
    text(zeros(height(sps),1),1:height(sps),num2str(round(sps.kept(:))),'horiz','center'); 
    text(-sps.lost(:)/2-sps.kept(:)/2,1:height(sps),num2str(round(sps.lost(:))),'horiz','center'); 
    text(sps.gain(:)/2+sps.kept(:)/2,1:height(sps),num2str(round(sps.gain(:))),'horiz','center'); 
    yline(0); box on; grid on; axis tight; xlim(x_lim)
end
% exportgraphics(gcf, "figures/phylo_ranking.eps")

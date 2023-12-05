

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

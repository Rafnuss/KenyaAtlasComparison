%% Load data
addpath("functions/")
sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");

name=["Marsh Warbler", "Sprosser", "Whitethroat", "River Warbler", "Irania", "Willow Warbler", "Red-backed Shrike", "Red-tailed Shrike", "Rufous Bush Chat", "Olive-tree Warbler", "Olivaceous Warbler", "Nightingale", "Basra Reed Warbler"];

[~,id] = ismember(name,sp.CommonName);

spj = sp(id,:);

spj.diff = (spj.new ./ spj.old).^(1/33);

figure;
bar(spj.diff)
xticks(1:height(spj))
xticklabels(spj.CommonName )
yline(1)
ylim([.95 1.05])

figure;box on; grid on;
b = barh([-spj.kept/2 -spj.lost spj.kept/2 spj.gain],1,'stacked',EdgeColor='none');
c=brewermap(3,'RdYlGn');
b(1).FaceColor = c(2,:);
b(2).FaceColor = c(1,:);
b(3).FaceColor = c(2,:);
b(4).FaceColor = c(3,:);
yticks(1:height(spj));
yticklabels(spj.CommonName);
% xlabel('Number of grid lost (-) and gainlost(+)')
set(gca, 'YDir','reverse')
text(zeros(height(spj),1),1:height(spj),num2str(round(spj.kept)),'horiz','center'); 
text(-spj.lost/2-spj.kept/2,1:height(spj),num2str(round(spj.lost)),'horiz','center'); 
text(spj.gain/2+spj.kept/2,1:height(spj),num2str(round(spj.gain)),'horiz','center'); 
grid on; axis tight; xlim([-90 90])


%% Turnover

tmp = (sp.lost+sp.gain)./(sp.lost+sp.gain+sp.kept);

figure; histogram(tmp)
hold on;
histogram(tmp(id))
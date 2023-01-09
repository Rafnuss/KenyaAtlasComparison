i_sp = find(strcmp(sp_old.CommonName,"Black Kite"));
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
im = imagesc(g.lon,g.lat,diff(:,:,i_sp),'alphadata',.8*(~isnan(diff(:,:,i_sp)))); 
tmp = diff(:,:,i_sp); tmp(tmp==0)=1;tmp=tmp.*C; id = ~isnan(tmp);
%scatter(g.LON(id),g.LAT(id),100,tmp(id)*8,'filled','MarkerEdgeColor','k')
%xline([g.lon(1)-g.res/2 g.lon+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
%yline([g.lat(1)-g.res/2 g.lat+g.res/2],Color=[.8 .8 .8], LineWidth=.5)
xticks([g.lon(1)-g.res/2 g.lon+g.res/2])
yticks([g.lat(1)-g.res/2 g.lat+g.res/2])
axis equal tight; set(gca,"YDir","normal");
% axis([33.4500   41.9500   -5.5500    5.0500])
title(sp_old.CommonName(i_sp))
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


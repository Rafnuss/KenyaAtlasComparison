%% Load data
addpath("functions/")
sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");
sp.diff = (sp.new - sp.old)./(sp.new + sp.old)*2;
sp = sp(sp.waterBird==1,:);

csr8 = readtable("data/CSR8/KE_slope_edit.csv", TextType="string");
csr8.meaning_sym(:) = "~";
csr8.meaning_sym(contains(csr8.meaning, "Strong decrease (p<0.01)")) = "⇊^{**}";
csr8.meaning_sym(contains(csr8.meaning, "Strong increase (p<0.01)")) = "⇈^{**}";
csr8.meaning_sym(contains(csr8.meaning, "Moderate decrease (p<0.01)")) = "↓^{**}";
csr8.meaning_sym(contains(csr8.meaning, "Moderate increase (p<0.01)")) = "↑^{**}";
csr8.meaning_sym(contains(csr8.meaning, "Strong decrease (p<0.05)")) = "⇊^*";
csr8.meaning_sym(contains(csr8.meaning, "Strong increase (p<0.05)")) = "⇈^*";
csr8.meaning_sym(contains(csr8.meaning, "Moderate decrease (p<0.05)")) = "↓^*";
csr8.meaning_sym(contains(csr8.meaning, "Moderate increase (p<0.05)")) = "↑^*";
csr8.meaning_sym(contains(csr8.meaning, "Stable")) = "→";


spj = outerjoin(sp, csr8, Keys="SEQ");


id = ~contains(spj.meaning, "Uncertain") & ~ismissing(spj.meaning);

%
fit(spj.diff(id), spj.add(id), fittype( @(a, x) a*x ), weight=(spj.new(id) + spj.old(id))/2)

corrW(spj.diff(id), spj.add(id), (spj.new(id) + spj.old(id))/2)
[rho,p] = corrcoef(spj.diff(id), spj.add(id));


figure; hold on; grid on;
scatter(spj.diff(id), spj.add(id),(spj.new(id) + spj.old(id))/2,'ok','filled')
errorbar(spj.diff(id), spj.add(id), spj.se_add(id),"xk","LineStyle","none")

l=lsline;
text(spj.diff(id), spj.add(id)+.001,spj.CommonName(id),"HorizontalAlignment","center")


%
figure; hold on; grid on;
tmp = spj.totals_imputed .* spj.add;
tmp = sign(tmp).*log(abs(tmp));
scatter(spj.new - spj.old, tmp, 10+10*log10(spj.totals_imputed),'ok','filled')
text(spj.new - spj.old, tmp,spj.CommonName,"HorizontalAlignment","center")


%
spj.meaning_sym(ismissing(spj.meaning_sym))="";
[~,tmp] = sort(spj.new - spj.old);
spj = spj(tmp,:);
figure('position',[0 0 3508 100+2480*height(spj)/100/2]/2); 
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
nd=66;
for i_s=1:3
    id_sub = (1:nd)+(i_s-1)*nd;
    id_sub = id_sub(id_sub<=height(spj));
    nexttile; box on; grid on;
    b = barh([-spj.kept(id_sub)/2 -spj.lost(id_sub) spj.kept(id_sub)/2 spj.gain(id_sub)],1,'stacked',EdgeColor='none');
    c=brewermap(3,'RdYlGn');
    b(1).FaceColor=  c(2,:);
    b(2).FaceColor=  c(1,:);
    b(3).FaceColor=  c(2,:);
    b(4).FaceColor=  c(3,:);
    yticks(1:height(spj));
    yticklabels(spj.CommonName(id_sub)+" "+spj.meaning_sym(id_sub));
    % xlabel('Number of grid lost (-) and gain(+)')
    set(gca, 'YDir','reverse')
    text(zeros(numel(id_sub),1),1:numel(id_sub),num2str(round(spj.kept(id_sub))),'horiz','center'); 
    text(-spj.lost(id_sub)/2-spj.kept(id_sub)/2,1:numel(id_sub),num2str(round(spj.lost(id_sub))),'horiz','center'); 
    text(spj.gain(id_sub)/2+spj.kept(id_sub)/2,1:numel(id_sub),num2str(round(spj.gain(id_sub))),'horiz','center'); 
    grid on; axis tight; xlim([-60 60])
end
exportgraphics(gcf, "export/traits/waterbirds.eps")
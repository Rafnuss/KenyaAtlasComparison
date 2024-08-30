
%% Load data
addpath("functions/")
sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");
% sp = sp((sp.lost+sp.gain+sp.kept)>10,:);
% sp(sp.CommonName=="House Sparrow",:) = [];
% sp(sp.CommonName=="Levant Sparrowhawk",:) = [];

burn = readtable("data/burns2021/EU birds decline overall in line with global patterns_species_results_withSEQ.csv", TextType="string");
burn.migratory_strategy = categorical(burn.migratory_strategy);
burn = burn(burn.SEQ>0,:);
spj = outerjoin(sp, burn, Keys="SEQ");


spj.diff = spj.gain_score - spj.lost_score;
spj.diff = (spj.new - spj.old)./(spj.new + spj.old)*2;
spj.diff = (spj.new ./ spj.old).^(1/33);
spj.diff = spj.new - spj.old;

id = spj.migratory_strategy~="Resident" & ~isnan(spj.annualRateOfChange) & ~isnan(spj.diff);% spj.migratory_strategy=="Resident"; spj.Palearctic==1 &
id = spj.palearctic==1 & ~isnan(spj.annualRateOfChange) & ~isnan(spj.diff) & (spj.lost+spj.gain+spj.kept)>10
id = ~isnan(spj.annualRateOfChange) & ~isnan(spj.diff) & (spj.lost+spj.gain+spj.kept)>1


fit(spj.diff(id), spj.annualRateOfChange(id),"poly1",weight=(spj.new(id) + spj.old(id))/2)
corrW(spj.diff(id), spj.annualRateOfChange(id), (spj.new(id) + spj.old(id))/2)
[rho,p] = corrcoef(spj.diff(id), spj.annualRateOfChange(id));

figure; hold on;
scatter(spj.diff(id), spj.annualRateOfChange(id),(spj.new(id) + spj.old(id)),'ok','filled')
l=lsline;
id2=id & spj.Migration~=3;% id2=id & spj.migratory_strategy~="Long-distance migrant";id2=id & spj.Palearctic~=1;
scatter(spj.diff(id2), spj.annualRateOfChange(id2),(spj.new(id2) + spj.old(id2)),'or','filled')

% text(spj.diff(id), spj.annualRateOfChange(id)+.001,spj.CommonName(id),"HorizontalAlignment","center")

yline(1)
xline(1)
box on; grid on; axis tight;
xlabel("Kenya - Average rate of change 1985 - 2010")
ylabel("Europe - Average rate of change 1980 - 2017")

%%
id2 = id & contains(spj.CommonName,"Warbler");
figure; hold on;
scatter(spj.diff(id2), spj.annualRateOfChange(id2),(spj.new(id2) + spj.old(id2)),'ok','filled')
text(spj.diff(id2), spj.annualRateOfChange(id2)+.001,spj.CommonName(id2),"HorizontalAlignment","center")

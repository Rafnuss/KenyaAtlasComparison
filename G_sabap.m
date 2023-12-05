
%% Compare raptor

addpath("functions/")
sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");
sp = sp((sp.lost+sp.gain+sp.kept)>10,:);
sp(sp.CommonName=="House Sparrow",:) = [];

sabap = readtable("data/SABAP1-2/underhill2014_withSEQ.csv", TextType="string");
sabap = sabap(sabap.SEQ>0,:);
sps = outerjoin(sp, sabap, Keys="SEQ");


z = [-1.96 -1.04 0 1.04 1.96];

tmp = sps{:,["z1","z2","z3","z4","z5","z6"]};
sps.diff_sabap = tmp * (1:6)' ./sum(tmp,2);


sps.diff = sps.new - sps.old;

figure; hold on,
scatter(sps.diff, sps.diff_sabap,sps.new + sps.old,'ok', 'filled')
lsline;
id = sps.Trophic_Niche=="Scavenger";
scatter(sps.diff(id), sps.diff_sabap(id),sps.new(id) + sps.old(id),'or', 'filled')
yticks(1:6); yticklabels(["less than -1.96", "-1.96 to -1.04", "-1.04 to 0", "0 to 1.04", "1.04 to 1.96", ">1.96"])

yline(3.5)
xline(0)
box on; grid on;
xlabel("Kenya - Change in Square 1985 - 2010")
ylabel("South Africa -  ﻿1987/91 - 2007/14")

corrW(sps.diff, sps.diff_sabap, (sps.new + sps.old)/2)

%%
id = sps.Trophic_Niche=="Scavenger";
sps(id,["CommonName_sabap", "diff_sabap"])
%% Compare raptor

addpath("functions/")
sp = readtable("export/sp_lost_kept_gain.csv", TextType="string");

D = table(repelem([99,91,89,92,87,129,88]',4,1),...
    repmat(["1980";"1980";"2000";"2000"],7,1),...    
    repmat(["non-migration";"migration"],7*2,1),...
    [3.04  6.84 4.1 7.78 0.4 7.23 0 0 24.14 116.76 21.27 46.1 2.17 7.21 2.30 1.26 4.96 14.16 4.44 5.23 6.19 8.3 5.93 4.43 1.83 1.87 1.36 0.7]',...
    VariableNames=["SEQ" "years" "season","v"]);

D_comp = unstack(D,"v","years");
D_comp.diff = (D_comp.x2000 - D_comp.x1980) ./ D_comp.x1980;
D_comp = unstack(removevars(D_comp,["x1980", "x2000"]),"diff","season");

[~,id] = ismember(D_comp.SEQ,sp.SEQ);

figure; hold on;
bar(categorical(sp.CommonName(id(id>0))), [D_comp.migration D_comp.non_migration])




sp_vult = sp(ismember(sp.SEQ,[D.SEQ; 90]),:);
% sp_vult.diff = (sp_vult.new ./ sp_vult.old).^(1/33);
sp_vult.diff = (sp_vult.new - sp_vult.old) ./ (sp_vult.new + sp_vult.old)*2;

sp_vult.diff = (sp_vult.new - sp_vult.old) ./ sp_vult.old;

sortrows(sp_vult(:,["CommonName","old","new","diff"]),"diff")



% sp(sp.Trophic_Niche=="Scavenger",:)



[]



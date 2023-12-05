

S = readtable("data/serengeti/dataset_modif.txt", TextType="string");
S.date = datetime(S.year,S.month,S.day);
S = removevars(S,["transectID" "year", "ID_" "month" "day" "specificEpithet" "macroHabitat" "observedDistance" "status" "family" "genus" "scientificNameAuthorship" "feedingLocation" "foodType" "locality" "observer"]);


sites = S(S.method=="Sites",:);

sites.site = string(sites.decimalLongitude) + "_" + string(sites.decimalLatitude);
sites.survey = sites.site+datestr(sites.date);

unique_survey = groupsummary(sites,["survey" "site" "date"]);
unique_survey.sum_species= unique_survey.GroupCount;
% figure; histogram(unique_survey.GroupCount)
unique_survey.before84 = year(unique_survey.date)<1985;

unique_sites = groupsummary(unique_survey,"site","sum",["sum_species" "before84"]);
unique_sites.sum_survey= unique_sites.GroupCount;

figure; histogram(unique_sites.sum_before84(unique_sites.sum_before84>0))

id = (unique_sites.sum_before84>=5) & (unique_sites.sum_survey-unique_sites.sum_before84)>5;



figure; hold on;
scatter(datenum(unique_survey.date),categorical(unique_survey.site),[],ismember(unique_survey.site, unique_sites.site(id)),'filled')
datetick("x")


%% 
% unique_sites.site(id)
i_site = "34.8491_-2.4379";

d = sites(sites.site==i_site,:);
d = removevars(d,["method" "decimalLongitude", "decimalLatitude" "site" "scientificName" "survey"]);

U = unstack(d,"individualCount","IOCEnglishName");
X = table2array(U(:,2:end))';

figure;
imagesc(X)
yticks(1:numel(U.Properties.VariableNames(2:end)))
yticklabels(U.Properties.VariableNames(2:end))
xticks(1:numel(U.date))
xticklabels(datestr(U.date))


figure; histogram(nansum(X),BinWidth=1)


any(~isnan(X(:,U.date<datetime(1984,1,1))),2)
any(~isnan(X(:,U.date>datetime(1984,1,1))),2)
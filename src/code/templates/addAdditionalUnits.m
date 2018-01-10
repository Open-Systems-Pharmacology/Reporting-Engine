
global MOBI_SETTINGS;
 
% make sure you application path is known and unit datie exist by initialisation of an exemplary xml file 
addpath(genpath('C:\Open Systems Pharmacology\Reporting Engine 1.0\'));
testReportingEngineInstallation('C:\Open Systems Pharmacology\Reporting Engine 1.0\');

% load unitFile
load(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'));

% check if desired Unit exist
jj = strcmp(unitList_dimensionList,'AUC (dose normalized)');
if ~any(jj)
    unitList_dimensionList{end+1} = 'AUC (dose normalized)';
    unitList(end+1) = unitList(end);

    unitList(end).unit_txt = {'h/l' 'min/l'};
    unitList(end).formula = {'60'  '1'};
    unitList(end).par_descriptions = {};
    unitList(end).par_names = {};
    unitList(end).par_values = [];
    unitList(end).baseUnit = 'min/l';
    
end

jj = strcmp(unitList_dimensionList,'AUC (BW normalized)');
if ~any(jj)
    unitList_dimensionList{end+1} = 'AUC (BW normalized)';
    unitList(end+1) = unitList(end);

    unitList(end).unit_txt = {'µg*h/L/kg'};
    unitList(end).formula = {'1'};
    unitList(end).par_descriptions = {};
    unitList(end).par_names = {};
    unitList(end).par_values = [];
    unitList(end).baseUnit = 'µg*h/L/kg';
    
end

jj = strcmp(unitList_dimensionList,'AUC (BSA normalized)');
if ~any(jj)
    unitList_dimensionList{end+1} = 'AUC (BSA normalized)';
    unitList(end+1) = unitList(end);

    unitList(end).unit_txt = {'µg*h/L/m²'};
    unitList(end).formula = {'1'};
    unitList(end).par_descriptions = {};
    unitList(end).par_names = {};
    unitList(end).par_values = [];
    unitList(end).baseUnit = 'µg*h/L/m²';
    
end

% sort alphabetical
[unitList_dimensionList,ix_sort]=sort(unitList_dimensionList);
unitList=unitList(ix_sort);

save(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'),'-append','unitList','unitList_dimensionList');

return


global MOBI_SETTINGS;
% initialize MoBiSettings
MoBiSettings;
 
% load unitFile
load(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'));

% check if desired Unit exist
iU = find(strcmp(unitList_dimensionList,'Amount'));
if ~any(strcmp(unitList(iU),'µmol'))
    % add uniq and formula for µg
    unitList(iU).unit_txt{end+1} = 'µmol';
    unitList(iU).formula{end+1} = '1/MW';
    % add uniq and formula for µg
    unitList(iU).unit_txt{end+1} = 'mg';
    unitList(iU).formula{end+1} = '1000/MW';
    % add pararmeter for translation to mass
    unitList(iU).par_descriptions = {'Molweight [g/mol]'};
    unitList(iU).par_names =  {'MW'};
    unitList(iU).par_values =  1;
end

iU = find(strcmp(unitList_dimensionList,'Mass'));
if ~any(strcmp(unitList(iU),'mg'))
    % add uniq and formula for µg
    unitList(iU).unit_txt{end+1} = 'mg';
    unitList(iU).formula{end+1} = 'MW';
    % add uniq and formula for µg
    unitList(iU).unit_txt{end+1} = 'µmol';
    unitList(iU).formula{end+1} = 'MW/1000';
    % add pararmeter for translation to mass
    unitList(iU).par_descriptions = {'Molweight [g/mol]'};
    unitList(iU).par_names =  {'µmol'};
    unitList(iU).par_values =  1;
end


save(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'),'-append','unitList','unitList_dimensionList');

return

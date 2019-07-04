function setAmountMassUnitList

global MOBI_SETTINGS;
% initialize MoBiSettings
if isempty(MOBI_SETTINGS)
    MoBiSettings;
end

% delete current unitFile
delete(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'));

[unitList,unitList_dimensionList]=iniUnitList(0);

% load clean unitFile
load(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'));

% Get Unit lists for Mass and Amount from Dimension
Mass_index = find(strcmp(unitList_dimensionList, 'Mass'));
Amount_index = find(strcmp(unitList_dimensionList, 'Amount'));

% Unify the unit list and formulas of amount and mass
Mass_unit_txt = [unitList(Mass_index).unit_txt, unitList(Amount_index).unit_txt];
% Unify the formulas of amount and mass
Mass_formula = [unitList(Mass_index).formula, strcat('1E-9*', unitList(Amount_index).formula, '*MW')];

% Unify the list of amount and mass
Amount_unit_txt = [unitList(Amount_index).unit_txt, unitList(Mass_index).unit_txt];
% Unify the formulas of amount and mass
Amount_formula = [unitList(Amount_index).formula, strcat('1E9*', unitList(Mass_index).formula, '/MW')];

%-----------------------------
% For adding Amount to Mass

% Check if the unit does not already contains the other values
if ~strContains(unitList(Mass_index).unit_txt, 'mol')
    % Update the values and formaulas for Mass
    unitList(Mass_index).par_descriptions = {'Molweight [g/mol]'};
    unitList(Mass_index).par_names = {'MW'};
    unitList(Mass_index).par_values = {'1'};
    
    unitList(Mass_index).unit_txt = Mass_unit_txt;
    unitList(Mass_index).formula = Mass_formula;
    
    save(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'),'-append','unitList','unitList_dimensionList');
end

%-----------------------------
% For add Mass to Amount

% Check if the unit does not already contains the other values
if ~strContains(unitList(Amount_index).unit_txt, 'kg')

    % Update the values and formaulas for Amount
    unitList(Amount_index).par_descriptions = {'Molweight [g/mol]'};
    unitList(Amount_index).par_names = {'MW'};
    unitList(Amount_index).par_values = {'1'};
    
    unitList(Amount_index).unit_txt = Amount_unit_txt;
    unitList(Amount_index).formula = Amount_formula;

    save(fullfile(MOBI_SETTINGS.application_path,'unitList_3.mat'),'-append','unitList','unitList_dimensionList');
end



return
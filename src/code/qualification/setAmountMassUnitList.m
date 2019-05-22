function unitList = setAmountMassUnitList(unitList, unitList_dimensionList)

% Get lists for Mass and Amount
Mass_index = find(strcmp(unitList_dimensionList, 'Mass'));
Amount_index = find(strcmp(unitList_dimensionList, 'Amount'));

% Unify the list of amount and mass
Mass_unit_txt = [unitList(Mass_index).unit_txt unitList(Amount_index).unit_txt];
Amount_unit_txt = [unitList(Amount_index).unit_txt unitList(Mass_index).unit_txt];

% Unify the formulas of amount and mass
Mass_formula = [unitList(Mass_index).formula strcat('1E-9*', unitList(Amount_index).formula, '*MW')];
Amount_formula = [unitList(Amount_index).formula strcat('1E9*', unitList(Mass_index).formula, '/MW')];

% Update the values and formaulas for Amount
unitList(Amount_index).par_descriptions = {'Molweight [g/mol]'};
unitList(Amount_index).par_names = {'MW'};
unitList(Amount_index).par_values = {'400'};

unitList(Amount_index).unit_txt = Amount_unit_txt;
unitList(Amount_index).formula = Amount_formula;

% Update the values and formaulas for Mass
unitList(Mass_index).par_descriptions = {'Molweight [g/mol]'};
unitList(Mass_index).par_names = {'MW'};
unitList(Mass_index).par_values = {'400'};

unitList(Mass_index).unit_txt = Mass_unit_txt;
unitList(Mass_index).formula = Mass_formula;


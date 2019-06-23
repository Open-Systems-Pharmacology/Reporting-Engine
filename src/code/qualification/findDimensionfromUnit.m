function [outputDimension, outputUnit, Flag] = findDimensionfromUnit(inputUnit)
% Find dimension from unit

% Load list of all units and dimensions
[unitList,unitList_dimensionList]=iniUnitList(-1, 3);

% Caution:
% Age and time share the same dimensions
% However, Time covers more unit and will be reference dimension
unitList(find(strcmp(unitList_dimensionList, 'Age in weeks'))).unit_txt=NaN;
unitList(find(strcmp(unitList_dimensionList, 'Age in years'))).unit_txt=NaN;

% Search for the 1st unit matching the input unit
for iUnit=1:length(unitList)
    
    MatchedUnit = find(strcmp(unitList(iUnit).unit_txt,inputUnit));
    
    if ~isempty(MatchedUnit)
        % Get the unit and dimension as cell
        outputUnit=inputUnit;
        outputDimension=unitList_dimensionList{iUnit};
        break
    end
end

Flag = isempty(MatchedUnit);

if Flag
    outputUnit=inputUnit;
    outputDimension=[];
end
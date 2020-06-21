function [outputDimension, outputUnit, Flag] = findDimensionfromUnit(inputUnit)
% Find dimension from unit

%fix the unit in case of wrong read or write encoding in cases where unit comes from file
inputUnit = adjustUnit(inputUnit);

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

function adjustedUnit = adjustUnit(unit)

try
    adjustedUnit = char(unit);
    
    %if the unit was e.g. read from file using wrong encoding, invalid characters are replace with characted code 65533
    %because in nearly 100% of cases invalid character is 'µ': Replace all occurances of char 65533 with 'µ'
    %If original character was not 'µ': exception will happen in any case (with or without replacement)
    adjustedUnit(double(adjustedUnit)==65533)='µ';

    %if the unit was e.g. read from file which was saved using wrong encoding char(194) is added before/after some characters. just remove it
    adjustedUnit = strrep(adjustedUnit,char(194),'');

catch
    %if something went wrong - return original unit
    adjustedUnit = unit;
end

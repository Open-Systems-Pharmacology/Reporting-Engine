function [outputPath, outputUnit, outputDimension] = strimPath(inputPath)
% Get unit and dimension from string path where units are within brackets
% Strim the path to get the output name

% Find in outputPath the unit within the brackets

OpeningBracketInd = strfind(inputPath,'[');
ClosingBracketInd = strfind(inputPath,']');

% If empty, there is no unit, output is dimensionless
if ~isempty(OpeningBracketInd)
    % In case more than one '[' ']' are found, takes the last
    OpeningBracketInd=OpeningBracketInd(end);
    ClosingBracketInd=ClosingBracketInd(end);
    
    % Split between header and unit
    outputPath = strtrim(inputPath(1:OpeningBracketInd-1));
    
    % Get unit within brackets
    tmpUnit=strtrim(inputPath(OpeningBracketInd+1:ClosingBracketInd-1)); %#ok<AGROW>
    
    try
        [outputDimension, outputUnit] = findDimensionfromUnit(tmpUnit);
    catch
        writeToReportLog('ERROR', sprintf('Unit %s extracted from path %s is not listed', tmpUnit, inputPath), 'false');
        error('Unit %s extracted from path %s is not listed', tmpUnit, inputPath);
    end
else
    outputPath=inputPath;
    outputUnit=[]; %#ok<AGROW>
    outputDimension='Fraction';
end
end
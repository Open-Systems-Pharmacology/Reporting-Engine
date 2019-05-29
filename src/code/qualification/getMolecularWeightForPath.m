function MW = getMolecularWeightForPath(outputPath,simulationIndex)
%GETMMOLECULARWEIGHTFORPATH Auxiliary function that get the molecular
% weight in g/l of a compound using path as input
%
% MW = getMolecularWeightForPath(outputPath,simulationIndex)
%
% Inputs:
%   outputPath (string) path to be parsed to get the compound
%   simulationIndex (integer), optional, index of the simulation
%
% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------
%

% check if simulationIndex is given

if ~exist('simulationIndex','var') || isempty(simulationIndex) || isnan(simulationIndex)
    
    simulationIndex = 1;
    
end


% get all Molecularweights and Compounds

[~,desc] = existsParameter('*|Molecular *',simulationIndex,'parameterType','readonly');

compounds = cellfun(@(x) x{end-1},regexp(desc(2:end,2),'\|','split'),'UniformOutput',false);

MW = cell2mat(desc(2:end,3));

% get compound which is part of path

jj = cellfun(@(x) contains([outputPath '|'],x),strcat('|',compounds,'|'));

% If getMolecularWeight does not find the compound in the path
% jj is a vector of zeros, then MW is NaN

if max(jj)==0
    
    MW = NaN;
    
else
    % If getMolecularWeight finds the compound in the path
    % set MW and output MW in g/l
    
    Compound = compounds{jj};
    MW = MW(jj);
    
    % Set MW in g/mol
    MWUnit = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly', 'property', 'Unit');
    
    % Issue with the encoding of molecular weight ? Output is 'kg/?mol' in
    % latest examples. Might be a version issue with the xml files ?
    try
        MW = MW.*getUnitFactor(MWUnit, 'g/mol', 'Molecular weight');
        
    catch
        [unitList,unitList_dimensionList]=iniUnitList(-1, 3);
        MWUnit = unitList(find(strcmp(unitList_dimensionList, 'Molecular weight'))).baseUnit;
        MW = MW.*getUnitFactor(MWUnit, 'g/mol', 'Molecular weight');
        
    end
end

return
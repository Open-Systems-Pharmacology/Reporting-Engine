function MW = getMolecularWeightForPath(outputPath,simulationIndex)

% check if simulationIndex is given

if ~exist('simulationIndex','var') || isempty(simulationIndex) || isnan(simulationIndex)
    
    simulationIndex = 1;
    
end

% get all Molecularweights

[~,desc] = existsParameter('*|Molecular *',simulationIndex,'parameterType','readonly');

compounds = cellfun(@(x) x{end-1},regexp(desc(2:end,2),'\|','split'),'UniformOutput',false);

MW = cell2mat(desc(2:end,3));

% get compound which is part of path

jj = cellfun(@(x) contains([outputPath '|'],x),strcat('|',compounds,'|'));

% set MW
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

return
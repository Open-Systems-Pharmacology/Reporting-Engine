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

%{
% If no or multiple compounds are found for that path,
% Boolean vector jj will lead to empty MW or multiple output
% In this case, the Molecular Weight entry is assumed
% and a warning output in the log

if sum(jj)==0
    jj=1;
    writeToReportLog('WARNING', sprintf(['In getMolecularWeightForPath: no match found between compounds and path %s \n' ...
        'Compound %s with MW= %f g/mol was assumed \n'], outputPath, compounds{jj}, MW(jj)), 'true');
elseif sum(jj)>1
    jj=find(jj==1, 1);
    writeToReportLog('WARNING', sprintf(['In getMolecularWeightForPath: multiple matches found between compounds and path %s \n' ...
        'Compound %s with MW= %f g/mol was assumed \n'], outputPath, compounds{jj}, MW(jj)), 'true');
end
%}

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
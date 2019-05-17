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

if sum(jj)~=1
    disp(jj);
    jj=1; 
end

% set MW
Compound = compounds{jj};
MW = MW(jj);

% Set MW in g/mol
MWUnit = getParameter(sprintf('*|%s|Molecular weight',Compound),1,'parametertype','readonly', 'property', 'Unit');

MW = MW.*getUnitFactor(MWUnit, 'g/mol', 'Molecular weight');


return
function [parValuesFinal,sourceIndex] = loadSelectedPathsOfPopulation(WSettings,listOfname,parPathSelection)
% LOADSELECTEDPATHSOFPOPULATION load values for selected paths for one or a merged population
% 
% parValuesFinal = loadSelectedPathsOfPopulation(WSettings,listOfname,parPathSelection)
% 
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       listOfname (cell array of strings)   list dfines the selection of  populations
%       parPathSelection (cell array of strings 2 x n) first row list the selected paths
%                               second row corresponding display units
% Outputs:
%      parValuesFinal (double matrix) (nInd x nPar) popParameter 
%      sourceIndex (double vector) (nInd x 1) index of source population

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


for iName = 1:length(listOfname)

    load(fullfile('tmp',listOfname{iName},'pop.mat'));
    [jj,ix] = ismember(parPathSelection(:,1),parPaths);
    
    if sum(jj) < size(parPathSelection,2)
       error('Error path selection for population parameter, was not correct');
    end

    if iName ==1
        % initialize parValuesFinal
        parValuesFinal = parValues(:,ix(jj)); %#ok<NODEF>
        sourceIndex = ones(size(parValuesFinal,1),1);
    else
        parValuesFinal = [parValuesFinal;parValues(:,ix(jj))]; %#ok<NODEF,AGROW>
        sourceIndex = [sourceIndex;ones(size(parValues,1),1).*iName]; %#ok<AGROW>
    end
end
        
% get Unit factor
unit = unit(ix(jj)); %#ok<NODEF>

for iU = 1:length(unit)
    
    switch unit{iU}
        case 'none'
            unitFactor = 1;
        otherwise
            unitFactor = getUnitFactorForUnknownDimension(WSettings,unit{iU},parPathSelection{iU,2});
    end
    
    parValuesFinal(:,iU) = parValuesFinal(:,iU).*unitFactor;

    
end
 

return

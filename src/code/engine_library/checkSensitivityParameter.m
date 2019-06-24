function sensParameterList = checkSensitivityParameter(sensParameterList,simulationIndex) 
%CHECKSENSITIVITYPARAMETER check parameter paths of sensitivity
%
% Inputs:
%       sensParameterList  (cellarry) 1. column pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. reportname
%                       5. column default value from xml
%                       6. minimal value
%                       7. maximal value
%       simulationIndex (double) indes of simulation


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% check if parameter name exist
keep = true(1,size(sensParameterList,1));
for iRow = 1:size(sensParameterList,1)
    if ~existsParameter(['*|' sensParameterList{iRow,1}],simulationIndex,'parametertype','readonly')
        writeToReportLog('WARNING',sprintf('%s is a parameter of the sensParameterList, which is not in the xml of the simulation %s, It will be omitted.',...
            sensParameterList{iRow}),false);
        keep(iRow) = false;
    else
        
        sensParameterList{iRow,5} = getParameter(['*|' sensParameterList{iRow,1}],simulationIndex,'parametertype','readonly');
        sensParameterList{iRow,6} = nan;
        sensParameterList{iRow,7} = nan;
    end
end
sensParameterList = sensParameterList(keep,:);

return



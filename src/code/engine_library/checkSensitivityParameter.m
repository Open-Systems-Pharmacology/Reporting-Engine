function sensParameterList = checkSensitivityParameter(WSettings,sensParameterList,simulationIndex)
%CHECKSENSITIVITYPARAMETER check parameter paths of sensitivity
%
% Inputs:
%       WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%       sensParameterList  (cellarry) 1. column pathid of parameter,
%                       2. number of steps
%                       3. variation range
%                       4. minimal value
%                       5. maximal value
%                       6. column default value from xml
%       simulationIndex (double) indes of simulation


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% check if parameter name exist
keep = true(1,size(sensParameterList,1));
for iRow = 1:size(sensParameterList,1)
    if ~existsParameter(['*|' sensParameterList{iRow,1}],simulationIndex,'parametertype','readonly')
        writeToLog(sprintf('WARNING %s is a parameter of the sensParameterList, which is not in the xml of the simulation %s, It will be omitted.',...
            sensParameterList{iRow}),WSettings.logfile,true,false);
        keep(iRow) = false;
    else
        
        sensParameterList{iRow,6} = getParameter(['*|' sensParameterList{iRow,1}],simulationIndex,'parametertype','readonly');
    end
end
sensParameterList = sensParameterList(keep,:);

return



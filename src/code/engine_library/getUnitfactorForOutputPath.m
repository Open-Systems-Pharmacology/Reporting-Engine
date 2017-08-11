function [unitFactor,success,MW] = getUnitfactorForOutputPath(WSettings,pathID,targetUnit,simulationIndex)
% GETUNITFACTORFOROUTPUTPATH getfactor to convert output defined by path to target unit
%
%  [unitFactor,success,MW] = getUnitfactorForOutputPath(WSettings,pathID,targetUnit,simulationIndex)
%
% Inputs: 
%      - WSettings (structure)    definition of properties used in all
%                   workflow functions see GETDEFAULTWORKFLOWSETTINGS
%      - pathID (string)  path which identifies output
%      - targetUnit (string) unit to convert internal unit
%      - simulationIndex (double) index of prevoiously initialized simulation
%  Outputs:
%      - unitFactor (double)  factor to convert interanl unit to target unit
%      - success (boolean) if false, it was not possible to determine unit conversion factor
%      - MW (double) molweigt of selected output

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



% initialie return values
unitFactor = nan;
success = true; %#ok<NASGU>
MW = nan; %#ok<NASGU>

% check if it is observer
ise =  existsObserver(['*' pathID],simulationIndex);
if ise
    unit = getObserverFormula(['*' pathID],simulationIndex,'property','Unit');

    % for concnetration the output should be normally an observer with the compound
    % name at second last path entry
    pathID = getObserverFormula(['*' pathID],simulationIndex,'property','Path');
    tmp = regexp(pathID,'\|','split');
    moleculeName = tmp{end-1};
    MW = getParameter(sprintf('*|%s|Molecular weight',moleculeName),simulationIndex,'parametertype','readonly');

else
     %check if it is statevariable
     ise = existsSpeciesInitialValue(['*' pathID],simulationIndex,'parametertype','readonly');
     if ise
         unit = getSpeciesInitialValue(['*' pathID],simulationIndex,'parametertype','readonly','property','Unit');
     end
     
     % species shhouldn't be in concentration
     MW = nan;
end
    
% if it neither observer nor statevariable, terminate function
if ~ise
    writeToLog(sprintf('ERROR: Outputpath "%s" could not be found in model',pathID),WSettings.logfile,true,false);
    success = false;
    return
end
    
[unitFactor,success] = getUnitFactorForUnknownDimension(WSettings,unit,targetUnit,MW);

return
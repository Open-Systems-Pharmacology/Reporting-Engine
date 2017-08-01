function [targetParameterValues]  = addDosetablePerBSA(WSettings,targetParameterList,parPaths,parValues,SplitCondition)
% ADDDOSETABLEPERBSA add new rows to the population, which overwrites the dose in the applicationprotocol
%
% [targetParameterValues]  = addDosetablePerBSA(WSettings,targetParameterList,parPaths,parValues,SplitCondition)
%
% Inputs:
%   - WSettings       (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   - targetParameterList (cellarray of string)   List of
%                       applicationparameter, which should be added whith adjusted dose to the
%                       population values
%   - parPaths (cellarray):  pathnames of the population parameters
%   - parValues (double matrix):  values of the population parameters and individuals
%   - SplitCondition structure with following fields
%           - BSAmin	 (double vector) [dm^2] vector which defines lower limit of one body weight bin
%           - BSAmax	 (double vector) [dm^2] vector which defines upper limit of one body weight bin
%           - targetParameter (double vector) values, which the target parameter will be in the corresponding body weight bin
%
% Outputs:
%   - targetParameterValues (double matrix):  values of the of the newly created population population parameters

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 18-July-2017


% get Body weight of poulation
jj = strcmp('Organism|Weight',parPaths);
weight =  parValues(:,jj);
jj = strcmp('Organism|Height',parPaths);
height =  parValues(:,jj);

BSA = calculateBodySurfaceArea(WSettings,weight,height);

% get new Values
targetParameterValues = nan(size(parValues,1), length(targetParameterList));

for iStep = 1:length(SplitCondition)
    
    jj = BSA >= SplitCondition(iStep).BSAmin & BSA < SplitCondition(iStep).BSAmax;
    targetParameterValues(jj,:) = SplitCondition(iStep).targetParameter;
    
end

% check if some individuals are outside bound
jj = BSA< min([SplitCondition.BSAmin]);
if any(jj)
    writeToLog(sprintf('%d individual(s) have/has a body surface area less than %g, and are not covered by the bsa dependent dosetable',...
        sum(jj),min([SplitCondition.BSAmin])),WSettings.logfile,true,false);
     targetParameterValues(jj,:) = 0;
end
jj = BSA >= max([SplitCondition.BSAmax]);
if any(jj)
    writeToLog(sprintf('%d individual(s) have/has a  body surface area greater equal than %g, and are not covered by the bsa dependent dosetable',...
        sum(jj),max([SplitCondition.BSAmax])),WSettings.logfile,true,false);

    targetParameterValues(jj,:) = 0;

end
    
return


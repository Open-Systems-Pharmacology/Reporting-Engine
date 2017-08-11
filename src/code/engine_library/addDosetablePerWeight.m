function [targetParameterValues]  = addDosetablePerWeight(WSettings,targetParameterList,parPaths,parValues,SplitCondition)
% add new lines to the population, which overwrites the dose in the applicationprotocol
%
% [targetParameterValues]  = addDosetablePerWeight(WSettings,targetParameterList,parPaths,parValues,SplitCondition)
%
% Inputs:
%   - WSettings       (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   - targetParameterList (cellarray of string)   List of
%                       applicationparameter, which should be added whith adjusted dose to the
%                       population values
%   - parPaths (cellarray):  pathnames of the population parameters
%   - parValues (double matrix):  values of the population parameters and individuals
%   - SplitCondition structure with following fields
%           - BWmin	 (double vector) [kg] vector which defines lower limit of one body weight bin
%           - BWmax	 (double vector) [kg] vector which defines upper limit of one body weight bin
%           - targetParameter (double vector) values, which the target parameter will be in the corresponding body weight bin
%
% Outputs:
%   - targetParameterValues (double matrix):  values of the of the newly created population population parameters

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% get Body weight of poulation
jjBW = strcmp('Organism|Weight',parPaths);
weight =  parValues(:,jjBW);


% get new Values
targetParameterValues = nan(size(parValues,1), length(targetParameterList));

for iStep = 1:length(SplitCondition)
    
    jj = weight >= SplitCondition(iStep).BWmin & weight < SplitCondition(iStep).BWmax;
    targetParameterValues(jj,:) = SplitCondition(iStep).targetParameter;
    
end

% check if some individuals are outside bound
jj = weight< min([SplitCondition.BWmin]);
if any(jj)
    writeToLog(sprintf('%d individual(s) have/has a bodyweight less than %g, and are not covered by the weight dependent dosetable',...
        sum(jj),min([SplitCondition.BWmin])),WSettings.logfile,true,false);
     targetParameterValues(jj,:) = 0;
end
jj = weight >= max([SplitCondition.BWmax]);
if any(jj)
    writeToLog(sprintf('%d individual(s) have/has a bodyweight greater equal than %g, and are not covered by the weight dependent dosetable',...
        sum(jj),max([SplitCondition.BWmax])),WSettings.logfile,true,false);

    targetParameterValues(jj,:) = 0;

end
    
return


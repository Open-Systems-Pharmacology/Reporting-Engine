function Settings = getDefaultWorkflowSettings
% GETDEFAULTWORKFLOWSETTINGS definition of properties used in all workflow functions
%  e.g. name of logfile, graphical properties, like colors which areused in all functions and can be customized
%
% Inputs - none
% Outputs Settings (structure)
% the fields of the structure contains two kinds of properties
% properties which can be customized
%   logfile     (string)    name of logfile, if empty it will be set by initialize workflow
%
% properties which are set internaly
%   isValidatedSystem   (boolean)   check if the workflow runs on a
%           validated system, otherwise all figures and tables are marked
%           with draft

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

%%  properties which can be customized
% name of logfile, if empty it will be set by initialize workflow to
% logfile.txt
Settings.logfile = '';

%% properties which are set internaly
% check if the workflow runs on a validated system
Settings.isValidatedSystem = false;


return




function [Settings] = initializeWorkflow(workflowType,Settings)
% initialise the workflow, sets up the logfile, and sets global settings
%
% Inputs:
%       workflowType (string)  type of workflow
%       Settings (structure)  see GETDEFAULTWORKFLOWSETTINGS
%
%

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 14-July-2017

%% adjust global settings
if isempty(Settings.logfile)
    Settings.logfile = 'logfile.txt';
end

% add Reporting engine version
TRE = ReportingEngineVersionInfo;


% getComputername
computerName = getenv('COMPUTERNAME');

% check if is a validated system 
Settings.isValidatedSystem = ismember(computerName,TRE.ListOfValidatedComputers);




% add OSPSuite version
T = OSPSuiteVersionInfo;
% T = SBSuiteVersionInfo;
if isunix %ToDo linux version
    versionOfSimMOdel = T.SimModelVersion; %default linux version
else
    versionOfSimMOdel = T.SimModelVersion_Linux; %default windows version
end

%% Initialize logfile
msg = sprintf('Start Workflow of type %s \n',workflowType);
if Settings.isValidatedSystem
    msg = sprintf('%s System is validated \n',msg);
else
    msg = sprintf('%s System is NOT validated \n',msg);
end
msg = sprintf('%s Computer %s (%s) \n',msg,computerName,computer);
msg = sprintf('%s Reporting Engine version: %g \n',msg,TRE.ReportingEngineVersion);
msg = sprintf('%s Open Systems Pharmacology Suite version: %g \n',msg,T.OSPSuiteVersion);
% msg = sprintf('%s Open Systems Pharmacology Suite version: %g \n',msg,T.SBSuiteVersion);
msg = sprintf('%s Version of SimModel: %g \n',msg,versionOfSimMOdel);
msg = sprintf('%s SimModelCompName: %s \n',msg,T.SimModelCompName);
msg = sprintf('%s Matlab version: %s \n',msg,version);

writeToLog(msg,Settings.logfile,true,true);

% close all precvious figures


return
function [WSettings] = initializeWorkflow(WSettings)
% initialise the workflow, sets up the logfile, and sets global settings
%
% [WSettings] = initializeWorkflow(workflowType,WSettings)
% 
% Inputs:
%       WSettings (structure)  see GETDEFAULTWORKFLOWSETTINGS
% Outputs:
%       WSettings (structure)  see GETDEFAULTWORKFLOWSETTINGS
%
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


%% adjust global settings

% add Reporting engine version
TRE = ReportingEngineVersionInfo;


% getComputername
computerName = getenv('COMPUTERNAME');

% check if is a validated system 
WSettings.isValidatedSystem = ismember(computerName,TRE.ListOfValidatedComputers);


% add OSPSuite version
T = OSPSuiteVersionInfo;
if isunix 
    versionOfSimMOdel = T.SimModelVersion_Linux; %default linux version
else
    versionOfSimMOdel = T.SimModelVersion; %default windows version
end

%% Initialize logfile
msg = sprintf('Start Workflow:  %s \n',workflowModeToText(WSettings.workflowType,WSettings.workflowMode));
if WSettings.isValidatedSystem
    msg = sprintf('%s System is validated \n',msg);
else
    msg = sprintf('%s System is NOT validated \n',msg);
end
msg = sprintf('%s Computer %s (%s) \n',msg,computerName,computer);
msg = sprintf('%s Reporting Engine version: %g \n',msg,TRE.ReportingEngineVersion);
msg = sprintf('%s Open Systems Pharmacology Suite version: %g \n',msg,T.OSPSuiteVersion);
msg = sprintf('%s Version of SimModel: %g \n',msg,versionOfSimMOdel);
msg = sprintf('%s SimModelCompName: %s \n',msg,T.SimModelCompName);
msg = sprintf('%s Matlab version: %s \n',msg,version);
msg = sprintf('%s Started by: %s \n',msg,getenv('Username'));

writeToReportLog('INFO',msg,true);
copyfile('logfile.txt','errorLogfile.txt');

% close all precvious figures
close all;

return
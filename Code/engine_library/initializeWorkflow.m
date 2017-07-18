function [Settings] = initializeWorkflow(workflowType,Settings)
% initialise the workflow, sets up the logfile, and sets global settings
%
% Inputs:
%       workflowType (string)  type of workflow
%       Settings (structure)  see GETDEFAULTWORKFLOWSETTINGS
%
%

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017


%% adjust global settings
if isempty(Settings.logfile)
    Settings.logfile = 'logfile.txt';
end

% check if is a vlaidated system % ToDo @Juri how can i do this customized
Settings.isValidatedSystem = false;


%% Initialize logfile
writeToLog(sprintf('Start Workflow of type %s \n computer %s \n Matlab version %s \n',...
    workflowType,computer,version),Settings.logfile,true,true); % %toDo @Juri how do i get the engine version?

return
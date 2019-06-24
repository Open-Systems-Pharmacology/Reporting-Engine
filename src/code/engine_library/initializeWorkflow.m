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

% if Watermark is not an input load default
if ~isfield(WSettings,'Watermark') 
    if WSettings.isValidatedSystem
        WSettings.Watermark = '';
    else
        userDir = fullfile(getenv('APPDATA') , SuiteInstallationSubfolder , 'ReportingEngine' );
        if ~exist(userDir,'dir')
            mkdir(userDir)
        end
        if exist(fullfile(userDir,'watermark.mat'),'file')
            load(fullfile(userDir,'watermark.mat')); %#ok<LOAD>
        else
           prompt={['Enter the default value for a watermark in plots' newline '(empty string will disable all watermarks):']};
           name='Watermark';
           numlines=1;
           defaultanswer={'Not QCed! Preliminary data!'};
 
           answer={};
           while isempty(answer)
              answer=inputdlg(prompt,name,numlines,defaultanswer);
           end
           watermark = answer{1};
           save(fullfile(userDir,'watermark.mat'),'watermark');
        end
        WSettings.Watermark = watermark;
    end
end
    


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

if ~WSettings.restart
    writeToReportLog('INFO',msg,true);
    copyfile('logfile.txt','errorLogfile.txt');
else
    writeToReportLog('INFO','Restart:',false);
    writeToReportLog('INFO',msg,false);
end
    

% close all precvious figures
close all;

return

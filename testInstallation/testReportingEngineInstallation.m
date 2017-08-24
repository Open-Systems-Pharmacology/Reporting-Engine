function isValid = testReportingEngineInstallation(rootFolderReportingEngine)

cd(fullfile(rootFolderReportingEngine,'testInstallation'))

% addpathes
addpath '..\src\code\engine_library';
addpath '..\src\matlab_toolbox\code';
addpath '..\src\matlab_toolbox\code\auxiliaries';

% name of xml
xml = 'Children_OralSingle_IV_Multi.xml';

% get defalt Setting
Settings = getDefaultWorkflowSettings('meanModel','default');

[~,isValid] = getApplicationProtocollFromXML(Settings,xml);

return
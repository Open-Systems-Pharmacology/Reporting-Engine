function isValid = testReportingEngineInstallation(rootFolderReportingEngine)

cd(fullfile(rootFolderReportingEngine,'testInstallation'))

% name of xml
xml = 'Children_OralSingle_IV_Multi.xml';

%init and process simulation
initSimulation(xml, []);
isValid=processSimulation;

return
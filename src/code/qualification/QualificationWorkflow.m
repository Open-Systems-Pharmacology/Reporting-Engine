% Script to start a Qualification Plan workflow
% Qualification Plan Workflow developed with Matlab 2017b
% --------------------------------------------------------------

clear
close all

% Set the working directy and create the output repository
% This part can be remove or replaced on later versions
REInput_path = 'QualificationPlan\examples\examples\advanced_01\reporting engine input';
cd(REInput_path)
REOutput_path = '..\reporting engine output';
mkdir(REOutput_path);

% --------------------------------------------------------------
% Get the Configuration Plan Settings
jsonFile = 'report-configuration-plan.json';
ConfigurationPlan = getConfigurationPlan(jsonFile);

% Linearize the Sections input and create the folder structure
ConfigurationPlan.Sections = generateOutputFolders(ConfigurationPlan.Sections, REOutput_path);

% --------------------------------------------------------------
% Copy Input content into right section
for i=1:length(ConfigurationPlan.Inputs)
    SectionPath = getSection(ConfigurationPlan.Sections, ConfigurationPlan.Inputs(i).SectionId);
    try
        copyfile(ConfigurationPlan.Inputs(i).Path, SectionPath)
        writeToReportLog('INFO',['Input ' num2str(i) ' was copied successfully'],'false');
    catch exception
        writeToReportLog('ERROR', ['Input ' num2str(i) ' could not be copied:' exception.message], 'true', exception);
        rethrow(exception);
    end
end


% Global settings
% Not developped yet
% WSettings = getDefaultWorkflowSettings('Qualification','default'); ??
WSettings = getDefaultWorkflowSettings('Qualification','default'); close;

% Load the Observations DataSets in structures
for i=1:length(ConfigurationPlan.ObservedDataSets)
    if ~iscell(ConfigurationPlan.ObservedDataSets)
        ConfigurationPlan.ObservedDataSets={ConfigurationPlan.ObservededDataSets};
    end
    
    ObservedDataSets(i) = loadObservationscsv(ConfigurationPlan.ObservedDataSets{1}.Path, ...
        ConfigurationPlan.ObservedDataSets{1}.Id);
end

% Task list: fields in ConfigurationPlan.Plots
TaskList = fields(ConfigurationPlan.Plots);

% runWorklfow
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

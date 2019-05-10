% Script to perform a Qualification Plan workflow
% Qualification Plan Workflow developed with Matlab 2017b
% --------------------------------------------------------------

clear
close all

% Set the working directy and create the output repository
% This part can be remove or replaced on later versions PKRatio advanced_01
% minimal PopulationTimeProfile
REInput_path = 'QualificationPlan\examples\minimal\reporting engine input';
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
% Copy Intro and Input content into right section
if isfield(ConfigurationPlan, 'Intro')
    for i=1:length(ConfigurationPlan.Intro)
        [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, ConfigurationPlan.Inputs(i).SectionId);
        try
            copyfile(ConfigurationPlan.Intro(i).Path, fullfile(REOutput_path));
        catch exception
            writeToReportLog('ERROR', sprintf('Input %d could not be copied: \n %s', num2str(i), exception.message), 'true', exception);
            rethrow(exception);
        end
    end
end
for i=1:length(ConfigurationPlan.Inputs)
    [SectionPath, indexed_item] = getSection(ConfigurationPlan.Sections, ConfigurationPlan.Inputs(i).SectionId);
    try
        copyfile(ConfigurationPlan.Inputs(i).Path, SectionPath);
    catch exception
        writeToReportLog('ERROR', sprintf('Input %d could not be copied: \n %s', num2str(i), exception.message), 'true', exception);
        rethrow(exception);
    end
end

% --------------------------------------------------------------
% Global settings
% Not developed yet: uses getDefaultWorkflowSettings from Reporting Engine
WSettings = getDefaultWorkflowSettings('Qualification','default'); close;

% --------------------------------------------------------------
% Load the Observations DataSets from csv into structures
if ~iscell(ConfigurationPlan.ObservedDataSets)
        ConfigurationPlan.ObservedDataSets=num2cell(ConfigurationPlan.ObservedDataSets);
end
for i=1:length(ConfigurationPlan.ObservedDataSets)
    % For DDI and PK Ratio plots, Observed Data Type is different
    % ObservedDataSets format to be determined
    % Temporary format: structure where the observed data is a table y
    if isfield(ConfigurationPlan.ObservedDataSets{i}, 'Type')
        ObsTable = readtable(ConfigurationPlan.ObservedDataSets{i}.Path, 'Encoding', 'UTF-8');
        ObservedDataSets(i).Id=ConfigurationPlan.ObservedDataSets{i}.Id;
        ObservedDataSets(i).Type=ConfigurationPlan.ObservedDataSets{i}.Type;
        ObservedDataSets(i).name=ConfigurationPlan.ObservedDataSets{i}.Path;
        ObservedDataSets(i).y=ObsTable;
    else
        ObservedDataSets(i) = loadObservationscsv(ConfigurationPlan.ObservedDataSets{i}.Path,...
            ConfigurationPlan.ObservedDataSets{i}.Id);
    end
end

% --------------------------------------------------------------
% Task list: fields in ConfigurationPlan.Plots
TaskList = fields(ConfigurationPlan.Plots);

% --------------------------------------------------------------
% run the Worklfow tasklist of ConfigurationPlan
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

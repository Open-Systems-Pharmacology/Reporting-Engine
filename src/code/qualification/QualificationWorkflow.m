% Script to start a Qualification Plan workflow
% Type: Qualification workflow
% Purpose:
% Original author: Pierre Chelle 01-01-2019
%
%  HOW TO USE
%
%  ! SPM reporting engine runs with Matlab 2013b (linux cluster)
%  Qualification plan developed with Matlab 2017b
% --------------------------------------------------------------

clear
close all

% Set the working directy
% This part can be remove or replaced on later versions
REInput_path = 'QualificationPlan\examples\minimal\reporting engine input';
cd(REInput_path)
REOutput_path = '..\reporting engine output';
mkdir(REOutput_path);

% Get the Configuration Plan Settings
jsonFile = 'report-configuration-plan.json';
ConfigurationPlan = getConfigurationPlan(jsonFile);

% Linearize the Sections input and create the folder structure
ConfigurationPlan.Sections = generateOutputFolders(ConfigurationPlan.Sections, REOutput_path);

% Global settings
% Not developped yet
% WSettings = getDefaultWorkflowSettings('Qualification','default'); ??
WSettings = getDefaultWorkflowSettings('',''); close;

% Load the Observations DataSets in structures
% Caution !! typo in name within json: ObservededDataSets
for i=1:length(ConfigurationPlan.ObservededDataSets)
    if ~iscell(ConfigurationPlan.ObservededDataSets)
        ConfigurationPlan.ObservedDataSets={ConfigurationPlan.ObservededDataSets};
    end
    
    ObservedDataSets(i) = loadObservationscsv(ConfigurationPlan.ObservededDataSets{1}.Path, ...
        ConfigurationPlan.ObservededDataSets{1}.Id);
end

% Task list: fields in ConfigurationPlan.Plots
TaskList = fields(ConfigurationPlan.Plots);

% runWorklfow
%runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,Datafiles,sensParameterList,MBS)
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

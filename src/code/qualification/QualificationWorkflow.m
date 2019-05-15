% Script to perform a Qualification Plan workflow
% Qualification Plan Workflow developed with Matlab 2017b
% --------------------------------------------------------------

clear all
close all

% List of examples to be tested 
REInput_path = 'QualificationPlanTests\examples\minimal\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\advanced_01\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\PKRatio\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\PopulationTimeProfile\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\advanced_01\reporting engine input';

REOutput_path = '..\reporting engine output';

cd(REInput_path)

% --------------------------------------------------------------
% Get the Configuration Plan Settings
jsonFile = 'report-configuration-plan.json';

[WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile, REOutput_path);

% --------------------------------------------------------------
% run the Worklfow tasklist of ConfigurationPlan
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

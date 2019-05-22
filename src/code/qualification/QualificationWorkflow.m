% Script to perform a Qualification Plan workflow
% Qualification Plan Workflow developed with Matlab 2017b
% --------------------------------------------------------------

clear all
close all
tic

% List of examples to be tested 
% REInput_path = 'QualificationPlanTests\examples\minimal\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\advanced_01\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\PKRatio\reporting engine input';
% REInput_path = 'QualificationPlanTests\examples\PopulationTimeProfile\reporting engine input';
<<<<<<< Updated upstream
REInput_path = 'C:\Design2Code\QualiExample01-master\re_input';
REOutput_path = 'C:\Design2Code\QualiExample01-master\re_output';
=======
<<<<<<< Updated upstream
REInput_path = 'QualificationPlanTests\examples\QualiExample01-master\re_input';
REOutput_path = 'QualificationPlanTests\examples\QualiExample01-master\re_output';
>>>>>>> Stashed changes

% REInput_path = 'QualificationPlanTests\examples\Qualification-Ontogeny-Distribution-GFR-master\re_input';

=======
 REInput_path = 'C:\Design2Code\QualiExample01-master\re_input';
 REOutput_path = 'C:\Design2Code\QualiExample01-master\re_output';
%REInput_path = 'C:\Design2Code\Qualification-Ontogeny-Distribution-GFR-master\re_input';
%REOutput_path = 'C:\Design2Code\Qualification-Ontogeny-Distribution-GFR-master\re_output';
>>>>>>> Stashed changes
% --------------------------------------------------------------
% Get the Configuration Plan Settings
jsonFile = 'report-configuration-plan.json';

[WSettings, ConfigurationPlan, TaskList, ObservedDataSets] = initializeQualificationWorkflow(jsonFile, REInput_path, REOutput_path);

% --------------------------------------------------------------
% run the Worklfow tasklist of ConfigurationPlan
runQualificationWorkflow(WSettings, ConfigurationPlan, TaskList, ObservedDataSets);

QualificationWorkflowTime = toc/60;
fprintf('\n Qualification Workflow Duration: %0.1f minutes \n', QualificationWorkflowTime);

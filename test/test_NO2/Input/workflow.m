% script to start the population workflow on the reporting engine
%
% DESCRIPTION
% a simulation is processed, PK Parameter are calculated for the selected
% timeprofiles and default plot for the simulations and population are
% generated
% simulated timeprofiles and PK Parmater are filed in the directory
% simulation, the format is suitable for an import to PK-sim
% the figures and tables in the directory figures
%
% HOW TO USE
% this script has to be filed in your working directory together with your input files like the simulation xml
% and the poulation csv
% adjust description for your purpose
% set the matlab director to your working directory
% start the script
%
% ! SPM reporting engine runs with Matlab 2013b (linux cluster)
%

% global settings 
% there are globale settings e.g. graphical properties, like colors which are
% used in all functions and can be customized
WSettings = getDefaultWorkflowSettings;

% generate workflowInput
[PopRunSet,TaskList,VPC] = generateWorkflowInputForPopulationSimulation(WSettings,'WorkflowInput.csv');


%% start the execution
runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC);
function [ TaskList ] = getDefaultTaskListPopulationWorkflow
% define tasklist for population workflow
%
% Inputs:    none
% Outputs:
%    - TaskList structure with following fields
%       simulatePopulation (boolean) if true the simulations defined in
%               popRunSet are simulated, false previously simulated results
%               are used (default true)
%       calculatePKParameter (boolean) if true the PK Parameters for simulations defined in
%               popRunSet are calculated, false previously simulated results
%               are used (default true)

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

% decides if the population has to be simulated or ob previously generated
% results should be used.
TaskList.simulatePopulation = true;

% decides if the PK-Parameter has to be calculated or ob previously generated
% results should be used.
TaskList.calculatePKParameter = true;


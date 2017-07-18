% ENGINE_LIBRARY
%
%  Master routines started by workflow scripts
%   runPopulationWorkflow                       - master routine for population workflow
%
%  functions called by master routines; which call one task specified
%   runPopulationSimulation                     - runs a population simulation defined within a PopRunSet 
%   calculatePopulationPKParameter              - calculates PK-Parameter for simulations defined within a PopRunSet 
%
%  functions, which defines default inputs used in the workflow script
%   getDefaultTaskListPopulationWorkflow        - define tasklist for population workflow
%   getDefaultPopRunSet                         - defines a  combinations of xml,population, outputs and study design to simulate
%   getDefaultWorkflowSettings                  - definition of properties used in all workflow functions
%
%  functions, which reads input files
%   readPopulationCSV                           - reads a poulation csv file 
%   readPopulationResultfile                    - reads a simulation result file exported by PK-SIM
%   readOutputSelection                         - reads list of selected outputs from csv file
%
%  internal functions
%   initializeWorkflow                          - initialise the workflow, sets up the logfile, and sets global settings
%   getApplicationProtocollFromXML              - get properties of the applicationsprotocol
%   calculatePKParameter_forApplicationProtocol - calculates PK Parameter, same List as in PK-Sim, for multi or single application

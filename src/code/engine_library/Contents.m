% ENGINE_LIBRARY
%
%  Master routines started by workflow scripts
%   runPopulationWorkflow                      - master routine for population workflow
%
%  functions called by master routines; which call one task specified
%   runPopulationSimulation                    - runs a population simulation defined within a PopRunSet 
%   calculatePopulationPKParameter             - calculates PK-Parameter for simulations defined within a PopRunSet 
%   runPopulationVPC                           - runs the visual predictive check for VPC
%
%  functions, which defines default inputs used in the workflow script
%   getDefaultTaskListPopulationWorkflow       - define tasklist for population workflow
%   getDefaultPopRunSet                        - defines a  combinations of xml,population, outputs and study design to simulate
%   getDefaultOutput                           - defines the properties of an Output
%   getDefaultWorkflowSettings                 - definition of properties used in all workflow functions
%   getDefaultVPCPopulationSettings            - get Settings for visual predicitve check
%
%  functions, which reads input files
%   readPopulationCSV                          - reads a poulation csv file 
%   readPopulationResultfile                   - reads a simulation result file exported by PK-SIM
%   readStudyDesign                            - read the study design description and adds this to the population
%   readNonmemFile                             - reads nonmemfile and executes filter,
%
%  functions, to create figures
%   plotReportHistogram                        - creates histogramm 
%
%  internal functions
%   initializeWorkflow                         - initialise the workflow, sets up the logfile, and sets global settings
%   getApplicationProtocollFromXML             - get properties of the applicationsprotocol
%   addDosetablePerWeight                      - add new lines to the population, which overwrites the dose in the applicationprotocol
%   calculatePKParameterForApplicationProtocol - calculates PK Parameter, same List as in PK-Sim, for multi or single application
%   getUnitfactorForOutputPath                 - getfactor to convert output defined by path to target unit
%   ReportFigurePrint                          - class which manages the print of figures
%   getReportFigure                            - creates new figure with and set watermark if not on validated system.
%   getUnitFactorForUnknownDimension           - getfactor to convert units, 
%   removeForbiddenLetters                     - replaces letters not suited for filenames with '_'
%   writeTabCellArray                          - write a cell array of strings as table
%
%  configuration function
%   ReportingEngineVersionInfo                 - configuration function, here properties of the reporting engine are listed
%  

% ENGINE_LIBRARY
%
%  Master routines started by workflow scripts
%   runPopulationWorkflow                        - master routine for population workflow
%
%  functions called by master routines; which call one task specified
%   runPopulationSimulation                      - runs a population simulation defined within a PopRunSet 
%   calculatePopulationPKParameter               - calculates PK-Parameter for simulations defined within a PopRunSet 
%   runPopulationVPC                             - runs the visual predictive check for VPC
%
%  functions, which defines default inputs used in the workflow script
%   getDefaultTaskListPopulationWorkflow         - define tasklist for population workflow
%   getDefaultOutput                             - defines the properties of an Output
%   getDefaultWorkflowSettings                   - definition of properties used in all workflow functions
%   getDefaultVPCPopulationSettings              - get WSettings for visual predicitve check
%
%  functions, which reads input files
%   readPopulationCSV                            - reads a poulation csv file 
%   readPopulationResultfile                     - reads a simulation result file exported by PK-SIM
%   readStudyDesign                              - read the study design description and adds this to the population
%   readNonmemFile                               - reads nonmemfile and executes filter,
%   readOutputCsv                                - reads the output csv and translate it to an array of structures
%   readPKParameter                              - reads a PK Analysis result file 
%
%  functions, to create figures
%   plotReportHistogram                          - creates histogramm 
%   plotReportShadedArea                         - Plots one property of a population vs another as shaed arae in comparison to a reference population
%   plotReportBoxwhisker                         - creates box whisker plot
%   plotReportTimeProfile                        - Plots the time profile of a population in comparison to a reference population
%
%  internal functions
%   generateWorkflowInputForPopulationSimulation - generates input according workflowinputcsv
%   initializeWorkflow                           - initialise the workflow, sets up the logfile, and sets global settings
%   getApplicationProtocollFromXML               - get properties of the applicationsprotocol
%   calculatePKParameterForApplicationProtocol   - calculates PK Parameter, same List as in PK-Sim, for multi or single application
%   plotLoopPopulationVPCPKparameter             - generates pK Parameter plots for VPC
%   plotLoopPopulationVPCforPhysiology           - does physiologocal plots for population VPC
%   plotLoopPopulationVPCTimeprofiles            - generates timeprofile plots for a population VPC
%
%  function which can be set by a function handle
%   addDosetablePerWeight                        - add new lines to the population, which overwrites the dose in the applicationprotocol
%   addDosetablePerBSA                           - add new rows to the population, which overwrites the dose in the applicationprotocol
%   textVPCPopulation                            - creates text fo figures and tables
%
%  functions to support graphical export
%   ReportFigurePrint                            - class which manages the print of figures
%   getReportFigure                              - creates new figure with and set watermark if not on validated system.
%
%  auxiliaries
%   getUnitfactorForOutputPath                   - getfactor to convert output defined by path to target unit
%   getUnitFactorForUnknownDimension             - getfactor to convert units, 
%   removeForbiddenLetters                       - replaces letters not suited for filenames with '_'
%   writeTabCellArray                            - write a cell array of strings as table
%   getLabelWithUnit                             - add unit to label if unit is not empty
%   calculateBodySurfaceArea                     - calculates the body surface area
%   loadSelectedPathsOfPopulation                - load values for selected paths for one or a merged population
%   loadMergedData                               - load data for different simualtion and merge them
%   getRangePlotPercentiles                      - calculates the upper mean and lower line for ranges plots
%   getcolmarkForMap                             - get colormatrix and marker-vector for a specified colormap
%   setLogarithmicYticks                         - set yticks to an axes, where the conten where plotted as log10
%
%  configuration function
%   ReportingEngineVersionInfo                   - configuration function, here properties of the reporting engine are listed
%  

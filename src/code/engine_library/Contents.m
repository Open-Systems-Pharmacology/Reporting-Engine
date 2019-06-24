% ENGINE_LIBRARY
%
%  Master routines started by workflow scripts
%   runPopulationWorkflow                      - master routine for population workflow
%   runMeanModelWorkflow                       - master routine for Mean Model workflow
%
%  functions called by master routines; which call one task specified
%   runPopulationSimulation                    - runs a population simulation defined within a PopRunSet 
%   calculatePopulationPKParameter             - calculates PK-Parameter for simulations defined within a PopRunSet 
%   runPopulationSensitivity                   - calculates and plots senstivity for all poulations
%   runPopulationVPC                           - runs the visual predictive check for VPC
%   runMeanModelVPC                            - runs the visual predictive check for VPC
%   runMeanModelAbsorption                     - gets absorption characteristic of the applications at time 0
%   runMeanModelSensitivity                    - runs the sensitivity analysis for all mean models defined in set
%   runMeanModelCheckMassbalance               - gets absorption characteristic of the applications at time 0
%
%  functions, which defines default inputs used in the workflow script
%   getDefaultWorkflowSettings                 - definition of properties used in all workflow functions
%   getDefaultVPCSettings                      - GETDEFAULTVPCPOPULATIONSETTINGS get Settings for visual predicitve check
%   getDefaultMassbalanceSettings              - GETDEFAULTVPCMEANMODELSETTINGS get WSettings for visual predicitve check
%
%
%  functions, which reads input files
%   readPopulationCSV                          - reads a poulation csv file 
%   readPopulationResultfile                   - reads a simulation result file exported by PK-SIM
%   readStudyDesign                            - read the study design description and adds this to the population
%   readNonmemFile                             - reads nonmemfile and executes filter,
%
%   readPKParameter                            - reads a PK Analysis result file 
%
%  functions, to create figures
%   plotReportHistogram                        - creates histogramm 
%   plotReportShadedArea                       - Plots one property of a population vs another as shaed arae in comparison to a reference population
%   plotReportBoxwhisker                       - creates box whisker plot
%   plotReportTimeProfile                      - Plots the time profile of a population in comparison to a reference population
%   plotReportPredictedVsObserved              - plot predcited vs observed
%   plotReportQQPlot                           - creates qqPlot for residuals
%   plotReportResiduals                        - plots residuals vs time or predicted y
%   plotSensListMostSensitive                  - plots all sensitivity above cutoff value
%
%  internal functions
%   initializeWorkflow                         - initialise the workflow, sets up the logfile, and sets global settings
%   getApplicationProtocollFromXML             - get properties of the applicationsprotocol
%   calculatePKParameterForApplicationProtocol - calculates PK Parameter, same List as in PK-Sim, for multi or single application
%   plotLoopPopulationVPCPKparameter           - generates pK Parameter plots for VPC
%   plotLoopPopulationVPCforPhysiology         - does physiologocal plots for population VPC
%   plotLoopVPCTimeprofiles                    - generates timeprofile plots for a population VPC
%   readTimeprofileDataForSimulation           - read timeprofile nonmen data and save the as temporary structure
%   calculatesPKParameterList                  - calculates the Pkparameter
%   generateSimResult                          - processes given  simulation with differnet parameter values
%   loadSimResult                              - load the simulated result for one sepcific output
%   calculateSensitivity                       - calculate sensitivity for PK Parameter
%   checkSensitivityParameter                  - check parameter paths of sensitivity
%   exportSensitivityForPopulation             - export list of all caclualted sensitivities
%   getListOfBestSensitivities                 - get all sensitivity above cutoff value
%   generateSensitivityParameterSet            - generates parameter for sensitivity calulation
%   addDependentPopulationParameter            - add additional population parameter 
%   getDefaultOutput                           - defines the properties of an Output
%
%
%  function which can be set by a function handle
%   addDosetablePerWeight                      - add new lines to the population, which overwrites the dose in the applicationprotocol
%   addDosetablePerBSA                         - add new rows to the population, which overwrites the dose in the applicationprotocol
%   textVPC                                    - TEXTVPCMEANMODEL creates text fo figures and tables
%
%
%  functions to support graphical export
%   ReportFigurePrint                          - class which manages the print of figures
%   getReportFigure                            - creates new figure with and set watermark if not on validated system.
%
%
%  auxiliaries
%   getUnitfactorForOutputPath                 - getfactor to convert output defined by path to target unit
%   getUnitFactorForUnknownDimension           - getfactor to convert units, 
%   removeForbiddenLetters                     - replaces letters not suited for filenames with '_'
%   writeTabCellArray                          - write a cell array of strings as table
%   getLabelWithUnit                           - add unit to label if unit is not empty
%   calculateBodySurfaceArea                   - calculates the body surface area
%   loadSelectedPathsOfPopulation              - load values for selected paths for one or a merged population
%   loadMergedData                             - load data for different simualtion and merge them
%   getRangePlotPercentiles                    - calculates the upper mean and lower line for ranges plots
%   getcolmarkForMap                           - get colormatrix and marker-vector for a specified colormap
%   setLogarithmicYticks                       - set yticks to an axes, where the conten where plotted as log10
%   checkInputDatafiles                        - check if varaibale Datafiles is correctly given
%   addToLegendPopulationSensitivity           - add legend entires to population sensitivity plots
%   writeToReportLog                           - Support function: Writes text to logfile. Each entry starts with a time stamp
%   workflowModeToText                         - generates text for logfile out of vartableworkflowMode
%   getVPCRange                                - calcualtes ranges for VPC plots
%   loadApplicationProtocoll                   - loads the temporary Outputlist
%   loadOutputList                             - loads the temporary Outputlist
%
%  configuration function
%   ReportingEngineVersionInfo                 - configuration function, here properties of the reporting engine are listed
%  

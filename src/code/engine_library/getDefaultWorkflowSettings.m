function WSettings = getDefaultWorkflowSettings(workflowType,workflowMode)
% GETDEFAULTWORKFLOWSETTINGS definition of properties used in all workflow functions
%  e.g. name of logfile, graphical properties, like colors which areused in all functions and can be customized
%
%  WSettings = getDefaultWorkflowSettings
% 
% Inputs 
%   workflowType (string) Type of workflow
%   
% Outputs 
%   WSettings (structure)
%       the fields of the structure contains two kinds of properties
%               properties which can be customized
%           logfile     (string)    name of logfile, if empty it will be set by initialize workflow
%           workflowType (string) type of workflow
%
%       properties which are set internaly
%           isValidatedSystem   (boolean)   check if the workflow runs on a
%                   validated system, otherwise all figures and tables are marked
%                   with draft

% todo documentation

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org



%%  properties which can be customized

% type of workflow
WSettings.workflowType = workflowType;
WSettings.workflowMode = workflowMode;

% caluclation Methodes
% Implemented are 'DuBois' and 'Haycock'
WSettings.BSA_calculationmethode = 'DuBois';

%% figure properties
% select output format possible
WSettings.printFormatList = {'EMF','FIG','CSV'};
% set name of figure directory
WSettings.figures = 'figures';

% define colors
%  (3x3) first color population, second color reference; third color data
colormap('default');
WSettings.colorVector = [[43,131,186]./255;  [171,221,164]./255;  [215,25,28]./255];
WSettings.colormapData = copper;
WSettings.colormapSimulation = winter;
WSettings.FaceAlpha = 0.5;

% used haded area and tiemprofile plots ( first and last) and for boxwhisker (all 5) 
WSettings.displayPercentiles = [5 25 50 75 95];

%Defines the methode for mean calclulation valid inputs are:
%       'median', 'geomean' or  'mean'  (arithmetic mean)
WSettings.shadedAreaMeanType = 'median';

% decide if boxwhsikers are plote with extremas
WSettings.boxwhiskerWithExtrema = false;
% angle of xlabel
WSettings.boxwhiskerXlabelRotation = 0;

% format for table export
WSettings.csvFormat = '%.3g';

% sensitivity analysis
WSettings.sensitivityCutoff = 0.9;
WSettings.sensitivitPrctileSelection = [5 50 95];

% number of indivduals per resultfile
WSettings.nIndPerSimResult = 1000;

% flag, if it is a restart should be used with care, if Matlab crashed
% during time consuming simulations, or a new run (default)
% restart, uses simulation and sensitivity results produced so far
WSettings.restart = false;

% minimal number for rangeplots
WSettings.rangePlotsMin = 2000;

% list of demographic parameters
WSettings.demographicParameters = {'IndividualId','Organism|Weight','Organism|Height','Organism|BMI'};

%% properties which are set internaly during inistialisation, so do not change them
% check if the workflow runs on a validated system
% getComputername
computerName = getenv('COMPUTERNAME');

% check if is a validated system % add Reporting engine version
TRE = ReportingEngineVersionInfo;
WSettings.isValidatedSystem = ismember(computerName,TRE.ListOfValidatedComputers);


return


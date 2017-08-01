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

% Open Systems Pharmacology Suite;  http://forum.open-systems-pharmacology.org
% Date: 28-July-2017

%% List all inputs files used for the workflow
% they should be filed in  your working directory
% for a population workflow you can generat the input files by exporting a
% populationsimulation "Export for Cluster Computations"
%
% list of all xml files, defining the simulations (source xxx.pksim5)
xmlList = {'SimplePop_Single_IVBolus.xml'};

% list of all population csv, 
% first column name of csv file, second column description used in figures and tables
popList = {'SimplePop_Single_IVBolus.csv','virtual population for test scripts'}; % (source xxx.pksim5)

% list of all outputs 
OutputList = {'Output.csv'};

% list of study design, like a dose table or application protocol, may be
% empty see exampleStudyDesign for more information
studyDesignList = {};


%% define the sets of population runs
% with a "PopRunSet"  you define which combinations of xml,population, outputs and
% study design you want to simulate,
% PopRunSet is a sructure see getDefaultPopRunSet
% inputs for getDefaultPopRunSet(name,reportName,xml,pop,OutputList)
% name (string) this name is used as directory or file name for the outputs, please avoid special signs like /
% reportName (string) this name is used in figures and tables
% xml (string) name of xml file
% pop (cell array of strings)    first entry name of csv file, second entry description for reporting
% OutputList (structure)  selection of prevoiusly defined outputs
PopRunSet(1)  =  getDefaultPopRunSet('SingleIvBolus','single IV application',xmlList{1},popList(1,:),OutputList{1}) ;

%% a population workflow can be started with different tasks
% speciy here which tasks should be executed in your case
TaskList = getDefaultTaskListPopulationWorkflow;

TaskList.simulatePopulation = false;
TaskList.calculatePKParameter = false;

TaskList.doVPC = true;
TaskList.doSensitivityAnalysis = false;

%% global settings 
% there are globale settings e.g. name of logfile, graphical properties, like colors which are
% used in all functions and can be customized
Settings = getDefaultWorkflowSettings;

%% get settings for Visualpredictive check
VPC = getDefaultVPCPopulationSettings(Settings,PopRunSet);

%% start the execution
runPopulationWorkflow(TaskList,Settings,PopRunSet,'VPC',VPC);
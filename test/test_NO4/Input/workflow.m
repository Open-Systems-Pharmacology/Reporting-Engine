% script which describes the pediatric workflow for a test example
% this script has to be filed in your working directory, and has to be
% started from the working directory, file all your other inputfiles also in
% this working directory, pathes will be set relative to this directory


%% List all inputs files used for the workflow
% they should be content of your working directory

% for a population workflow you can generat the input files by exporting a
% populationsimulation "Export for Cluster Computations"
% list of all xml files, defineing the simulations
xmlList = {'Children_OralSingle_IV_Multi.xml'};
% list of all population csv, 
% first column name of csv file, second column description for reporting
popList = {'Children_OralSingle_IV_Multi.csv','virtual population for test scripts'};

% list of all outputs 
outputcsv = {'Output.csv'};

% list of study design, like a dose table or application protocol
studyDesignList = {};

% name of nonmemfile containing the timeprofile data, 
%set to '', if no data available
dataTpFile = '';

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
PopRunSet(1)  =  getDefaultPopRunSet('OralSingle_IV_Multi','multi IV and PO application',xmlList{1},popList(1,:),outputcsv{1}) ;


%% a population workflow can be started with different tasks
% speciy here which tasks should be executed in your case
TaskList = getDefaultTaskListPopulationWorkflow;

TaskList.doVPC = false;
TaskList.doSensitivityAnalysis = false;

%% global settings 
% there are glbale settings e.g. fieldnames selected out of nonmemfile or graphical properties, like colors which are
% used in all functions and can be customized
Settings = getDefaultWorkflowSettings;


%% start the execution
runPopulationWorkflow(TaskList,Settings,PopRunSet);
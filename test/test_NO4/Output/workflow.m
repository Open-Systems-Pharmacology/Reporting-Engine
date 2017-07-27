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

% list of all outputs defined as structure
OutputList(1) = getDefaultOutput('Organism|PeripheralVenousBlood|Hydroxy_Itraconazole|Plasma (Peripheral Venous Blood)',...
    'Hydroxy Itraconazole','µg/l');
OutputList(2) = getDefaultOutput('Organism|PeripheralVenousBlood|Itraconazole|Plasma (Peripheral Venous Blood)',...
    'Itraconazole','µg/l');
OutputList(3) = getDefaultOutput('Organism|PeripheralVenousBlood|Midazolam|Plasma (Peripheral Venous Blood)',...
    'Midazolam','µg/l');
OutputList(4) = getDefaultOutput('Organism|Bone|Intracellular|Midazolam|Concentration in container',...
    'Midazolam Bone','µg/l');
OutputList(5) = getDefaultOutput('Organism|Bone|Intracellular|Itraconazole|Concentration in container',...
    'Itraconazole Bone','µg/l');
OutputList(6) = getDefaultOutput('Organism|Bone|Intracellular|Hydroxy_Itraconazole|Concentration in container',...
    'Hydroxy Itraconazole Bone','µg/l');


% list of study design, like a dose table or application protocol
studyDesignList = {};

% name of nonmemfile containing the timeprofile data, 
%set to '', if no data available
dataTpFile = '';

%% define the sets of population runs
% with a set you define which combinations of xml,population, outputs and
% study design you want to simulate, 
% this name is used as directory or file name for the outputs, please avoid special signs like /
simulationName = 'OralSingle_IV_Multi';
PopRunSet(1)  =  getDefaultPopRunSet(simulationName,xmlList{1},popList(1,:),OutputList(1:6)) ;


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
% Script to start a workflow
% Type: Population workflow, with project specific mode parallelComparison
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 23-Aug-2017 18:16:48
% 
%  HOW TO USE
%  this script has to be filed in your working directory together with your input files like the simulation xml
%  and the poulation csv
%  adjust description for your purpose
%  set the matlab directory to your working directory
%  start the script
% 
%  ! SPM reporting engine runs with Matlab 2013b (linux cluster)
% 


% global settings
% there are globale settings which are used in all functions.
WSettings = getDefaultWorkflowSettings('popModel','parallelComparison');

% Definitions of sets of Populations
clear PopRunSet
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|Theophylline|Plasma (Peripheral Venous Blood)','reportName','Theophylline','displayUnit','mg/l','dataTpFilter','CMT==1','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC','Thalf','CL','Vss';'mg/l','h','mg/l','mg*h/l','h','ml/min/kg','ml/kg'};
PopRunSet(1) = struct('name','PO320mg','reportName','PO administration of 320 mg','boxwhiskerLabel','Theophylline','xml','Theophylline (Boeckmann 1994) PO for fitting (fitted specific int.perm) popsim.xml','studyDesign','','isReference',0,'popcsv','Theophylline (Boeckmann 1994) PO for fitting (fitted specific int.perm) popsim.csv','popReportName','virtual population','calculatePKParameterFh','','dataTpFilter','STUD==9999','dataReportName','Study data 9999','OutputList',OutputList);
clear OutputList

% Definitions of TaskList
TaskList = struct('simulatePopulation',1,'calculatePKParameter',1,'doVPC',1,'doSensitivityAnalysis',0);

% List of Nonmemfiles:
% set to {}, if no data available
%  first column nonmem file, second column dictionary, third column datatype
dataFiles(1,:) = {'Theophylline_Nmdat.mndat','tpDictionary.csv','timeprofile'};

% get Definition of default plots
if TaskList.doVPC
    VPC = getDefaultVPCSettings(WSettings,PopRunSet);
else
    VPC = [];
end

% List of Parameters for sensitivity analysis:
% set to {}, if not needed
%  columns: 1. path, 2. number of steps, 3. variation range, 4. minValue 5. maxValue
sensParameterList = {};

% start the execution
runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,dataFiles,sensParameterList);

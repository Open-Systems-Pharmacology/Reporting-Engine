% Script to start a workflow
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 10-Aug-2017 16:20:43
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
WSettings = getDefaultWorkflowSettings;

% Definitions of sets of Populations
clear PopRunSet
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC','CL','Vss';'µg/l','h','µg/l','µg*h/l','ml/min/kg','ml/kg'};
OutputList(2) = struct('pathID','Organism|Bone|Intracellular|C1|Concentration in container','reportName','C1 bone intracellular','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
PopRunSet(1) = struct('name','SingleIvBolus','reportName','single IV application','boxwhiskerLabel','','xml','SimplePop_Single_IVBolus.xml','studyDesign','','isReference',0,'popcsv','SimplePop_Single_IVBolus.csv','popReportName','virtual test population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList

% Definitions of TaskList
TaskList = struct('simulatePopulation',1,'calculatePKParameter',1,'doVPC',1,'doSensitivityAnalysis',1);

% List of Nonmemfiles:
% set to {}, if no data available
%  first column nonmem file, second column dictionary, third column datatype
dataFiles = {};

% get Definition of default plots
if TaskList.doVPC
    VPC = getDefaultVPCPopulationSettings(PopRunSet,'parallelComparison');
else
    VPC = [];
end

% start the execution
runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,dataFiles);

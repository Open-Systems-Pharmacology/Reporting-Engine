% Script to start a workflow
% Type: Population workflow, with project specific mode default
% Purpose:: Test Case 3 Test PK Parameter calculation with function handle
% M&S activity:: None
% Validation level:: None
% Original author: ZTCOK 14-Nov-2017 13:00:32
% 
%  HOW TO USE
%  this script has to be filed in your working directory together with your input files like the simulation xml
%  and the poulation csv
%  adjust description for your purpose
%  set the matlab directory to your working directory
%  start the script
% 
%  Script is intendended to run with Matlab 2017b 
% 


% global settings
% there are globale settings which are used in all functions.
WSettings = getDefaultWorkflowSettings('popModel','default');

% Definitions of sets of Populations
clear PopRunSet
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','C_max_norm','t_max','C_tEnd','AUC','AUC_norm','AUC_inf','AUC_inf_norm','MRT','Thalf','FractionAucLastToInf','CL','Vss','Vd';'µmol/l','µg/l','h','µg/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','','ml/min/kg','ml/kg','ml/kg'};
OutputList(2) = struct('pathID','Organism|Bone|Intracellular|C1|Concentration in container','reportName','C1 bone intracellular','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(2).pKParameterList = {'C_max','C_max_norm','t_max','C_tEnd','AUC','AUC_norm','AUC_inf','AUC_inf_norm','MRT','Thalf','FractionAucLastToInf','CL','Vss','Vd';'µg/l','mg/l','h','µg/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','','ml/min/kg','ml/kg','ml/kg'};
PopRunSet(1) = struct('name','SingleIvBolus','reportName','single IV application','boxwhiskerLabel','','xml','SimplePop_Single_IVBolus.xml','studyDesign','','isReference',0,'popcsv','SimplePop_Single_IVBolus.csv','popReportName','virtual test population','calculatePKParameterFh','myCalculatePKParameterForApplicationProtocol','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList

% Definitions of TaskList
TaskList = struct('simulatePopulation',0,'calculatePKParameter',1,'doVPC',0,'doSensitivityAnalysis',0);

% List of Nonmemfiles:
% set to {}, if no data available
%  first column nonmem file, second column dictionary, third column datatype
dataFiles = {};

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

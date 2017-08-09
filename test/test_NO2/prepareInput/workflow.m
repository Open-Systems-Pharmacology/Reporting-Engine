% Script to start a workflow
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 09-Aug-2017 16:13:23
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
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|Hydroxy_Itraconazole|Plasma (Peripheral Venous Blood) with Typo','reportName','Hydroxy Itraconazole','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','C_max_norm','C_max_t1_t2','C_max_t1_t2_norm','C_max_tLast_tEnd','C_max_tLast_tEnd_norm','t_max','t_max_t1_t2','t_max_tLast_tEnd','C_trough_t2','C_trough_tLast','AUC_t1_t2','AUC_t1_t2_norm','AUC_tLast_minus_1_tLast','AUC_tLast_minus_1_tLast_norm','AUC_inf_t1','AUC_inf_t1_norm','AUC_inf_tLast','AUC_inf_tLast_norm','MRT','Thalf','Thalf_tLast_tEnd';'µmol/l','mg/l','µmol/l','mg/l','µmol/l','mg/l','h','h','h','µmol/l','µmol/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','h'};
OutputList(2) = struct('pathID','Organism|PeripheralVenousBlood|Itraconazole|Plasma (Peripheral Venous Blood)','reportName','Itraconazole','displayUnit','cm','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(2).pKParameterList = {'C_max','C_max_norm','C_max_t1_t2','C_max_t1_t2_norm','C_max_tLast_tEnd','C_max_tLast_tEnd_norm','t_max','t_max_t1_t2','t_max_tLast_tEnd','C_trough_t2','C_trough_tLast','AUC_t1_t2','AUC_t1_t2_norm','AUC_tLast_minus_1_tLast','AUC_tLast_minus_1_tLast_norm','AUC_inf_t1','AUC_inf_t1_norm','AUC_inf_tLast','AUC_inf_tLast_norm','MRT','Thalf','Thalf_tLast_tEnd';'µmol/l','mg/l','µmol/l','mg/l','µmol/l','mg/l','h','h','h','µmol/l','µmol/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','h'};
OutputList(3) = struct('pathID','Organism|PeripheralVenousBlood|Midazolam|Plasma (Peripheral Venous Blood)','reportName','Midazolam','displayUnit','typo','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(3).pKParameterList = {'C_max','C_max_norm','C_max_t1_t2','C_max_t1_t2_norm','C_max_tLast_tEnd','C_max_tLast_tEnd_norm','t_max','t_max_t1_t2','t_max_tLast_tEnd','C_trough_t2','C_trough_tLast','AUC_t1_t2','AUC_t1_t2_norm','AUC_tLast_minus_1_tLast','AUC_tLast_minus_1_tLast_norm','AUC_inf_t1','AUC_inf_t1_norm','AUC_inf_tLast','AUC_inf_tLast_norm','MRT','Thalf','Thalf_tLast_tEnd';'µmol/l','mg/l','µmol/l','mg/l','µmol/l','mg/l','h','h','h','µmol/l','µmol/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','h'};
OutputList(4) = struct('pathID','Organism|Bone|Intracellular|Itraconazole|Concentration in container','reportName','Itraconazole Bone','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(4).pKParameterList = {'C_max','C_max_norm','C_max_t1_t2','C_max_t1_t2_norm','C_max_tLast_tEnd','C_max_tLast_tEnd_norm','t_max','t_max_t1_t2','t_max_tLast_tEnd','C_trough_t2','C_trough_tLast','AUC_t1_t2','AUC_t1_t2_norm','AUC_tLast_minus_1_tLast','AUC_tLast_minus_1_tLast_norm','AUC_inf_t1','AUC_inf_t1_norm','AUC_inf_tLast','AUC_inf_tLast_norm','MRT','Thalf','Thalf_tLast_tEnd';'µmol/l','typo2','µmol/l','mg/l','µmol/l','mg/l','h','h','h','µmol/l','µmol/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','h'};
OutputList(5) = struct('pathID','Organism|Bone|Intracellular|Hydroxy_Itraconazole|Concentration in container','reportName','Hydroxy Itraconazole Bone','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(5).pKParameterList = {'C_max','C_max_norm','C_max_t1_t2','C_max_t1_t2_norm','C_max_tLast_tEnd','C_max_tLast_tEnd_norm','t_max','t_max_t1_t2','t_max_tLast_tEnd','C_trough_t2','C_trough_tLast','AUC_t1_t2','AUC_t1_t2_norm','AUC_tLast_minus_1_tLast','AUC_tLast_minus_1_tLast_norm','AUC_inf_t1','AUC_inf_t1_norm','AUC_inf_tLast','AUC_inf_tLast_norm','MRT','Thalf','Thalf_tLast_tEnd';'h','mg/l','µmol/l','mg/l','µmol/l','mg/l','h','h','h','µmol/l','µmol/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','µmol*min/l','µg*min/l','h','h','h'};
PopRunSet(1) = struct('name','OralSingle_IV_Multi','reportName','Multi IV and PO application','boxwhiskerLabel','','xml','Children_OralSingle_IV_Multi.xml','studyDesign','','isReference',0,'popcsv','Children_OralSingle_IV_Multi.csv','popReportName','virtual test population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
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

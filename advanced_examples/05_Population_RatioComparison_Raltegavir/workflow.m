% Script to start a workflow
% Type: Population workflow, mode ratioComparison
% Purpose:: Building and evaluation of a PBPK model for raltegravir in children
% M&S activity:: Raltegravir
% Validation level:: B
% Original author: XXXXX 07-Feb-2020 13:24:15
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
WSettings = getDefaultWorkflowSettings('popModel','ratioComparison');

% Definitions of sets of Populations
clear PopRunSet
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|Raltegravir|Plasma (Peripheral Venous Blood)','reportName','Raltegravir','displayUnit','µg/l','dataTpFilter','strcmp(OUTPUT,"Raltegravir_PLASMA")','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','C_max_norm','t_max','C_tEnd','AUC','AUC_norm','AUC_inf','AUC_inf_norm','MRT','Thalf','CL','Vss','Vd';'µg/l','kg/l','h','µg/l','µg*h/l','kg*h/l','µg*h/l','kg*h/l','h','h','ml/min/kg','ml/kg','ml/kg'};
PopRunSet(1) = struct('name','Sick_adults_400mg_sd','reportName','Sick adults 400mg sd','boxwhiskerLabel','Sick adults 400mg sd','xml','Raltegravir 400mg filmcoated tablet.xml','studyDesign','','isReference',0,'popcsv','Raltegravir Sick Adult Population.csv','popReportName','Sick Adult Population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|Raltegravir|Plasma (Peripheral Venous Blood)','reportName','Raltegravir','displayUnit','µg/l','dataTpFilter','strcmp(OUTPUT,"Raltegravir_PLASMA")','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','C_max_norm','t_max','C_tEnd','AUC','AUC_norm','AUC_inf','AUC_inf_norm','MRT','Thalf','CL','Vss','Vd';'µg/l','kg/l','h','µg/l','µg*h/l','kg*h/l','µg*h/l','kg*h/l','h','h','ml/min/kg','ml/kg','ml/kg'};
PopRunSet(2) = struct('name','Healthy_adults_400mg_sd','reportName','Healthy adults 400mg sd','boxwhiskerLabel','Healthy adults 400mg sd','xml','Raltegravir 400mg filmcoated tablet.xml','studyDesign','','isReference',1,'popcsv','Raltegravir Adult Population.csv','popReportName','Healthy Adult Population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList

% Definitions of TaskList
TaskList = struct('simulatePopulation',1,'calculatePKParameter',1,'doVPC',1,'doSensitivityAnalysis',0);

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
sensParameterList(1,:) = {'Raltegravir|Lipophilicity',2,0.1,'Lipophilicity'};
sensParameterList(2,:) = {'Raltegravir|Specific intestinal permeability (transcellular)',2,0.1,'Specific intestinal permeability (transcellular)'};
sensParameterList(3,:) = {'Raltegravir-UGT1A9-Kassahun 2007|In vitro Vmax for liver microsomes',2,0.1,'Raltegravir-UGT1A9-Kassahun 2007|In vitro Vmax for liver microsomes'};
sensParameterList(4,:) = {'Raltegravir-UGT1A1-Kassahun 2007|In vitro Vmax for liver microsomes',2,0.1,'Raltegravir-UGT1A1-Kassahun 2007|In vitro Vmax for liver microsomes'};
sensParameterList(5,:) = {'Applications|*|filmcoated tablet (original Merck formulation)|Dissolution shape',2,0.1,'Dissolution shape'};

% start the execution
runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,dataFiles,sensParameterList);

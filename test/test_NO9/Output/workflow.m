% Script to start a workflow
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 09-Aug-2017 17:03:02
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
PopRunSet(1) = struct('name','adult','reportName','single IV administration','boxwhiskerLabel','adult','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Adults_18_100_years-Population.csv','popReportName','virtual adult population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC','CL','Vss';'µg/l','h','µg/l','µg*h/l','ml/min/kg','ml/kg'};
PopRunSet(2) = struct('name','children_0_3y','reportName','single IV administration','boxwhiskerLabel','0-3y','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Kids_0_3_years-Population.csv','popReportName','virtual children 0-3 years','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC','CL','Vss';'µg/l','h','µg/l','µg*h/l','ml/min/kg','ml/kg'};
PopRunSet(3) = struct('name','children_3_6y','reportName','single IV administration','boxwhiskerLabel','3-6y','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Kids_3_6_years-Population.csv','popReportName','virtual children 3-6 years','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC','CL','Vss';'µg/l','h','µg/l','µg*h/l','ml/min/kg','ml/kg'};
PopRunSet(4) = struct('name','children_6_9y','reportName','single IV administration','boxwhiskerLabel','6-9y','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Kids_6_9_years-Population.csv','popReportName','virtual children 6-9 years','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
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

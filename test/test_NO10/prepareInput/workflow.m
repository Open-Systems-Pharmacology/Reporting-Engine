% Script to start a workflow
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 16-Aug-2017 15:10:04
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
OutputList(1).pKParameterList = {'C_max','AUC';'µg/l','µg*h/l'};
PopRunSet(1) = struct('name','adult','reportName','single IV administration','boxwhiskerLabel','adult','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',1,'popcsv','Adults_18_100_years-Population.csv','popReportName','virtual adult population','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','AUC';'µg/l','µg*h/l'};
PopRunSet(2) = struct('name','children_0_3y','reportName','single IV administration','boxwhiskerLabel','0-3y','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Kids_0_3_years-Population.csv','popReportName','virtual children 0-3 years','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','AUC';'µg/l','µg*h/l'};
PopRunSet(3) = struct('name','children_3_6y','reportName','single IV administration','boxwhiskerLabel','3-6y','xml','Adults_SingleIV_Bolus.xml','studyDesign','','isReference',0,'popcsv','Kids_3_6_years-Population.csv','popReportName','virtual children 3-6 years','calculatePKParameterFh','','dataTpFilter','','dataReportName','','OutputList',OutputList);
clear OutputList
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|C1|Plasma (Peripheral Venous Blood)','reportName','C1 plasma','displayUnit','µg/l','dataTpFilter','','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','AUC';'µg/l','µg*h/l'};
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
    VPC = getDefaultVPCPopulationSettings(PopRunSet,'pediatric');
else
    VPC = [];
end

% List of Parameters for sensitivity analysis:
% set to {}, if not needed
%  columns: 1. path, 2. number of steps, 3. variation range, 4. minValue 5. maxValue
sensParameterList(1,:) = {'Organism|Hematocrit',2,0.1,0,1};
sensParameterList(2,:) = {'Organism|Bone|Volume',2,0.1,0,NaN};
sensParameterList(3,:) = {'Organism|Fat|Volume',2,0.1,0,NaN};
sensParameterList(4,:) = {'Organism|Fat|Fraction interstitial',2,0.1,0,1};
sensParameterList(5,:) = {'Organism|Stomach|Volume',2,0.1,0,NaN};
sensParameterList(6,:) = {'Organism|Stomach|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(7,:) = {'Organism|SmallIntestine|Volume',2,0.1,0,NaN};
sensParameterList(8,:) = {'Organism|SmallIntestine|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(9,:) = {'Organism|LargeIntestine|Volume',2,0.1,0,NaN};
sensParameterList(10,:) = {'Organism|LargeIntestine|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(11,:) = {'Organism|Liver|Fraction vascular',2,0.1,0,1};
sensParameterList(12,:) = {'Organism|Liver|Fraction interstitial',2,0.1,0,1};
sensParameterList(13,:) = {'Organism|Liver|Volume',2,0.1,0,NaN};
sensParameterList(14,:) = {'Organism|Liver|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(15,:) = {'Organism|Pancreas|Volume',2,0.1,0,NaN};
sensParameterList(16,:) = {'Organism|Pancreas|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(17,:) = {'Organism|Spleen|Volume',2,0.1,0,NaN};
sensParameterList(18,:) = {'Organism|Spleen|Specific blood flow rate',2,0.1,0,NaN};
sensParameterList(19,:) = {'Neighborhoods|Bone_int_Bone_cell|BAY 1007626|P_cell_int_factor',2,0.1,0,NaN};
sensParameterList(20,:) = {'Neighborhoods|Bone_int_Bone_cell|BAY 1007626|P_int_cell_factor',2,0.1,0,NaN};
sensParameterList(21,:) = {'Neighborhoods|Fat_int_Fat_cell|BAY 1007626|P_cell_int_factor',2,0.1,0,NaN};

% start the execution
runPopulationWorkflow(WSettings,TaskList,PopRunSet,'pediatric',VPC,dataFiles,sensParameterList);

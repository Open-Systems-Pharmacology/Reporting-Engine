% Script to start a workflow
% Purpose:
% M&S activity:
% Validation level:
% Original author: ZTCOK 10-Aug-2017 16:19:21
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

% Definitions of sets of Mean Model silulations
clear MeanModelSet
OutputList(1) = struct('pathID','Organism|PeripheralVenousBlood|Theophylline|Plasma (Peripheral Venous Blood)','reportName','Theophylline','displayUnit','mg/l','dataTpFilter','STUD==0','residualScale','log','unitFactor',NaN,'pKParameterList',[]);
OutputList(1).pKParameterList = {'C_max','t_max','C_tEnd','AUC';'mg/l','h','mg/l','mg*h/l'};
MeanModelSet(1) = struct('name','PO320mg','reportName','PO administration of 320 mg','xml','Theophylline (Boeckmann 1994) PO.xml','calculatePKParameterFh','','dataTpFilter','STUD==0','dataReportName','Boeckmann 1994','OutputList',OutputList);
clear OutputList

% Definitions of TaskList
TaskList = struct('doVPC',1,'doSensitivityAnalysis',0);

% List of Nonmemfiles:
% set to {}, if no data available
%  first column nonmem file, second column dictionary, third column datatype
dataFiles = {'data.mndat','dict.csv','timeprofile'};

% get Definition of default plots
if TaskList.doVPC
    VPC = getDefaultVPCMeanModelSettings(MeanModelSet);
else
    VPC = [];
end

% start the execution
runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,dataFiles);

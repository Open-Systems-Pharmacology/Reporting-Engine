function [PopRunSet,TaskList,VPC] = generateWorkflowInputForPopulationSimulation(WSettings,workflowInputCsv)
% GENERATEWORKFLOWINPUTFORPOPULATIONSIMULATION generates input according workflowinputcsv
%
% Inputs 
%   WSettings           (structure) containing global settings see GETDEFAULTWORKFLOWSETTINGS
%   workflowInputCsv    (string) Path of inout csv, which defines the workflow 
% 
% Outputs
%    - PopRunSet        (structure) with following fields
%       - name  (string) as input variable name
%       - xml   (string) as input variable xml
%       - OutputList    (struture)  List of outputs see GETDEFAULTOUTPUT
%       - studyDesign   (string) csv which defines study design like dose
%                   table this is an optional input, it may be empty
%       - dataTpFilter  (string) this is an optional input, it may be empty
%               if timeprofile data are available, and a filter per popRunSet is necessary 
%               please enter here the filter conditions. 
%               use nonmem column headers as variable name.
%               e.g. 'AGE>4 & STUD==11111'
%       - isReference   (boolean)  marks if it is a refernce population 
%                   (e.g. adult in pediatric or healthy in disease comparison, 
%                   which shall be treated differently in plotting routines
%       - calculatePKParameterFH (function handle) 
%               if empty defualt function is used
%
%    - TaskList structure with following fields
%       simulatePopulation (boolean) if true the simulations defined in
%               popRunSet are simulated, false previously simulated results
%               are used 
%       calculatePKParameter (boolean) if true the PK Parameters for simulations defined in
%               popRunSet are calculated, false previously simulated results
%               are used 
%       doVPC (boolean) if true the VPC plots are generated 
%       doSensitivityAnalysis (boolean) if true a sensitivity analysis is done  


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% read data
data = readtab(workflowInputCsv,';',0,0,1,0);


%% generate Population sets:

% get rows with identifiers for simulation
jjIdentifier = ~cellfun(@isempty,data{3,1}) & ~cellfun(@(x) strncmp(x,'Task',4),data{3,1}) ...
    & ~cellfun(@(x) strncmp(x,'Workflow',8),data{3,1});

for iD = 3:size(data,2)
    
    % get values of simulation
    tmp = data{3,iD}(jjIdentifier);
    
    % convert booleans
    jj0 = cellfun(@(x) strcmp(x,'0'),tmp);
    tmp(jj0) = {false};
    jj1 = cellfun(@(x) strcmp(x,'1'),tmp);
    tmp(jj1) = {true};
    
    PopRunSet(iD-2) = cell2struct(tmp,data{3,1}(jjIdentifier)); %#ok<AGROW>
end

%% generate Tasklist

% get rows with identifiers for Tasks
jjIdentifier = ~cellfun(@isempty,data{3,1}) & cellfun(@(x) strncmp(x,'Task',4),data{3,1}) ...
    & ~cellfun(@(x) strncmp(x,'Workflow',8),data{3,1});

% get values of first simulation
tmp = data{3,3}(jjIdentifier);

% convert booleans
jj0 = cellfun(@(x) strcmp(x,'0'),tmp);
tmp(jj0) = {false};
jj1 = cellfun(@(x) strcmp(x,'1'),tmp);
tmp(jj1) = {true};

% get header
header = strrep(data{3,1}(jjIdentifier),'Task','');
% convert to struct
TaskList = cell2struct(tmp,header);

%% get Workflowtype
% get rows with identifiers for Tasks
jjIdentifier = strcmp('WorkflowType',data{3,1});
workflowType = data{3,3}{jjIdentifier};

if ~ismember(workflowType,{'pediatric','parallelComparison'})
    error('workflowType %s is invalid',workflowType);
end


% get VPC
if TaskList.doVPC
    VPC = getDefaultVPCPopulationSettings(WSettings,PopRunSet,workflowType);
else
    VPC = [];
end
function [SimulationSet,TaskList,Workflow,dataFiles,sensParameterList,outputSheets] = readWorkflowInput(workflowInputxls,sheet)
% READWORKFLOWINPUT converts workflow xls to STructures
%
%  Inputs 
%       workflowInputxls (string) name of xls file
%       sheet (string) name of sheet, if empty first sheet is taken
%  Outputs 
%       SimulationSet (structure) stores property of MeanModelset or RunPoulationSet
%       TaskList (structure) stores list of tasks
%       workflowType (string)  type of workflow
%       dataFiles (cellarray) strores name of data files
%       sensParameterList (cellarray) definiton for sensitivity analysis

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% read data
[~,~,xlsdata] = xlsread(workflowInputxls,sheet);

xlsIdentifier = xlsdata(2:end,1);
xlsInfo = xlsdata(2:end,3:end);
% remove empty rows
jj = ~cellfun(@isnumeric,xlsIdentifier) & ~cellfun(@isempty,xlsIdentifier);
xlsIdentifier = xlsIdentifier(jj);
xlsInfo = xlsInfo(jj,:);

% select and sort rows
% rows containing a referenced sheet
jjSheet =  strncmp(xlsIdentifier,'sheet',5);
% rows containing a Task
jjTask =  strncmp(xlsIdentifier,'Task',4);
% rows containing a Workflowdescription
jjWorkflow = strncmp(xlsIdentifier,'Workflow',8);
% rows with identifiers for simulation
jjIdentifier =  ~jjTask  & ~jjSheet & ~jjWorkflow...
    & ~ismember(xlsIdentifier,{'dataFileTimeprofile'});

%% check if refrences sheets are available
% get rows wirh sheet references

tmp = xlsInfo(jjSheet,:);
tmp = reshape(tmp,size(tmp,1)*size(tmp,2),1);
jj = ~cellfun(@isnumeric,tmp);
if any(jj)
    referencedSheets = unique(tmp(jj));

    [~,xlsSheets] = xlsfinfo(workflowInputxls);

    jj = ~ismember(referencedSheets,xlsSheets);
    if any(jj)
        error('Missing sheets: %s',strjoin(referencedSheets(jj),'; '))        
    end
end
%% generate Simulation Set:
% get cols for simulation info
jjCol =  ~cellfun(@isempty,xlsInfo(strcmp(xlsIdentifier,'name'),:)) & ...
     ~cellfun(@isnumeric,xlsInfo(strcmp(xlsIdentifier,'name'),:));


iD = 0;
for iCol = find(jjCol)
    
    iD = iD +1;
    
    % get values of simulation
    tmp = xlsInfo(jjIdentifier,iCol);

    % set nan to empty
    ij = find(cellfun(@isnumeric,tmp));
    jj = cellfun(@isnan,tmp(ij));
    tmp(ij(jj)) = {''};

    % convert booleans
    jj0 = cellfun(@(x) strcmp(x,'0'),tmp);
    tmp(jj0) = {false};
    jj1 = cellfun(@(x) strcmp(x,'1'),tmp);
    tmp(jj1) = {true};
    
    SimulationSet(iD) = cell2struct(tmp,xlsIdentifier(jjIdentifier)); %#ok<AGROW>
end
% check for unique identifier and replace special signs
if length(unique({SimulationSet.name})) < length(SimulationSet)
    error( 'Name of your simulations are not unique')
end
forbiddenLetters = {' ','.','/','\'};
for iD = 1:length(SimulationSet)
    for iL = 1:length(forbiddenLetters)
        SimulationSet(iD).name = strrep(SimulationSet(iD).name,forbiddenLetters{iL},'_');
    end
end
        


%% read outputs
outputSheets = xlsInfo(strcmp('sheetOutput',xlsIdentifier),jjCol);

%% generate Tasklist
% get values of first simulation
tmp = xlsInfo(jjTask,1);

if any(cellfun(@isempty,tmp))
    error('Please inset 0 or 1 to Tasks')
end
% convert booleans
jj0 = cellfun(@(x) strcmp(x,'0'),tmp);
tmp(jj0) = {false};
jj1 = cellfun(@(x) strcmp(x,'1'),tmp);
tmp(jj1) = {true};

% delete word "Task" to create header
header = strrep(xlsIdentifier(jjTask),'Task','');
% convert to struct
TaskList = cell2struct(tmp,header);


%% get Workflow Description

tmp = xlsInfo(jjWorkflow,1);
% delete word "Workflow" to create header
header = strrep(xlsIdentifier(jjWorkflow),'Workflow','');
% convert to struct
Workflow = cell2struct(tmp,header);

if ~isfield(Workflow,'Mode')
    Workflow.Mode = 'default';
else
    if ~ismember(Workflow.Mode,{'pediatric','parallelComparison','ratioComparison'})
        error('workflowMode %s is invalid',Workflow.Mode);
    end
end

%% generate Datafiles
[~,ijIdentifier] = ismember({'dataFileTimeprofile','sheetDataDictTimeprofile'},xlsIdentifier);
dataFiles = xlsInfo(ijIdentifier,1);
% all fields are filled
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
if ~all(cellfun(@isnumeric,dataFiles))
    dataFiles = {dataFiles{1}, [dataFiles{2} '.csv'],'timeprofile'};
    t = readtable(workflowInputxls,'Sheet',dataFiles{2}(1:end-4));
    writetable(t,dataFiles{2},'Delimiter',';','WriteRowNames',false);
    % only one is filled
elseif any(cellfun(@isnumeric,tmp))
    error('If you want to use data, fill all data fields, otherwise leave them empty')
else
    dataFiles = {};
end
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');

%% get sensitivity
sensSheet = xlsInfo{strcmp('sheetSensitivity',xlsIdentifier),1};
if isnumeric(sensSheet)
    sensParameterList={};
else
    [~,~,sensParameterList] = xlsread(workflowInputxls,sensSheet);
    if size(sensParameterList,2)<4
        error('Please insert 4 columns for sensitivities: Path, number of steps, variation range,report name');
    end
    sensParameterList = sensParameterList(2:end,1:4);
    jj = ~cellfun(@isnumeric,sensParameterList(:,1));
    sensParameterList = sensParameterList(jj,:);
end

return

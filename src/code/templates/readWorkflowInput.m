function [SimulationSet,TaskList,workflowType,dataFiles,sensParameterList] = readWorkflowInput(workflowInputxls,sheet)
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


%% generate Population sets:

% get rows with identifiers for simulation
jjIdentifier =  ~strncmp(xlsIdentifier,'Task',4) ...
    & ~strncmp(xlsIdentifier,'sens',4) ...
    & ~ismember(xlsIdentifier,{'WorkflowType','dataFileTimeprofile','dataDictTimeprofile'});

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


% check if outputxls is empty
if isfield(SimulationSet,'outputxls')
    jj = cellfun(@isempty,{SimulationSet.outputxls}) & ...
        ~cellfun(@isempty,{SimulationSet.outputsheet});
    
    if any(jj)
        for iSet  = find(jj)
            SimulationSet(iSet).outputxls = workflowInputxls; %#ok<AGROW>
        end
    end
end





%% generate Tasklist

% get rows with identifiers for Tasks
jjIdentifier = ~cellfun(@isempty,xlsIdentifier) & strncmp(xlsIdentifier,'Task',4);


% get values of first simulation
tmp = xlsInfo(jjIdentifier,1);

% convert booleans
jj0 = cellfun(@(x) strcmp(x,'0'),tmp);
tmp(jj0) = {false};
jj1 = cellfun(@(x) strcmp(x,'1'),tmp);
tmp(jj1) = {true};

% get header
header = strrep(xlsIdentifier(jjIdentifier),'Task','');
% convert to struct
TaskList = cell2struct(tmp,header);


%% get Workflowtype
workflowType = '';
% get rows with identifiers for Tasks
jjIdentifier = strcmp('WorkflowType',xlsIdentifier);
if any(jjIdentifier)
    workflowType = xlsInfo{jjIdentifier,1};
    
    if ~ismember(workflowType,{'pediatric','parallelComparison','ratioComparison'})
        error('workflowType %s is invalid',workflowType);
    end
end

%% generate Datafiles
dataFiles = {};
[jjIdentifier,ijIdentifier] = ismember({'dataFileTimeprofile','dataDictTimeprofile'},xlsIdentifier);
if all(jjIdentifier)
    tmp = xlsInfo(ijIdentifier,1);
    if ~all(cellfun(@isnumeric,tmp));
        dataFiles = [tmp', {'timeprofile'}] ;
    end
end

%% get sensitivity
[jjIdentifier,ijIdentifier] = ismember({'sensXls','sensSheet'},xlsIdentifier);
if all(jjIdentifier)
    tmp = xlsInfo(ijIdentifier,1);
    jj = cellfun(@isnumeric,tmp);
    if all(jj)  &&  all(cellfun(@isnan,tmp(jj)))
        sensParameterList={};
    else
        if isnan(tmp{1})
            sensXLS = workflowInputxls;
        else
            sensXLS = tmp{1};
        end
        if isnan(tmp{2})
            sensSheet = 1;
        else
            sensSheet = tmp{2};
        end
        [~,~,sensParameterList] = xlsread(sensXLS,sensSheet);
        if size(sensParameterList,2)<4
            error('Please insert 4 columns for sensitivities: Path, number of steps, variation range,report name');
        end
        sensParameterList = sensParameterList(2:end,1:4);
        jj = ~cellfun(@isnumeric,sensParameterList(:,1));
        sensParameterList = sensParameterList(jj,:);
    end
end

return

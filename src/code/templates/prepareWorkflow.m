function prepareWorkflow(workflowInputxls)
% PREPAREPOPULATIONWORKFLOW creates out of xls table workflow script
%
%  Inputs 
%   workflowInputxls (string) name of xls file
%   sheet (string) name of sheet, if empty first sheet is taken

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

% read Documentation
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
documentation = readtable(workflowInputxls,'Sheet','Documentation','ReadVariableNames',false);
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');

[SimulationSet,TaskList,Workflow,dataFiles,sensParameterList,outputSheets] = readWorkflowInput(workflowInputxls,1);


%% add Outputs
for iSet = 1:length(SimulationSet)
    SimulationSet(iSet).OutputList = readOutputXls(workflowInputxls,outputSheets{iSet});
end


%%write new worflow.m
fid = fopen('workflow.m','w');
% write header
% write header
writeWorkflowHeader(fid,Workflow.Type,Workflow.Mode,documentation);

switch Workflow.Type
    case 'popModel'
        writePopulations(fid,SimulationSet);
    case 'meanModel'
        writeMeanModel(fid,SimulationSet);
    otherwise
        error('unknown workflowtype')
end


% write Definitions of sets of TaskList
writeTaskList(fid,TaskList);


% write Datafiles
writeDataFiles(fid,dataFiles);

%% write function call VPC
fprintf(fid,'%% get Definition of default plots');
fprintf(fid,'\r\n');
fprintf(fid,'if TaskList.doVPC');
fprintf(fid,'\r\n');
switch Workflow.Type
    case 'popModel'
        fprintf(fid,'    VPC = getDefaultVPCSettings(WSettings,PopRunSet);');
    case 'meanModel'
        fprintf(fid,'    VPC = getDefaultVPCSettings(WSettings,MeanModelSet);');       
end
fprintf(fid,'\r\n');
fprintf(fid,'else');
fprintf(fid,'\r\n');
fprintf(fid,'    VPC = [];');
fprintf(fid,'\r\n');
fprintf(fid,'end');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

switch Workflow.Type
    case 'popModel'
   
    case 'meanModel'
        %% write function call MBS
        fprintf(fid,'%% get Definition of default Masbalance plots');
        fprintf(fid,'\r\n');
        fprintf(fid,'if TaskList.checkMassbalance');
        fprintf(fid,'\r\n');
        fprintf(fid,sprintf('    MBS = getDefaultMassbalanceSettings;'));
        fprintf(fid,'\r\n');
        fprintf(fid,'else');
        fprintf(fid,'\r\n');
        fprintf(fid,'    MBS = [];');
        fprintf(fid,'\r\n');
        fprintf(fid,'end');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\n');
end

%% write function call Sensitivity
% write sensitivity Aalysis
writeSensitivityParameterList(fid,sensParameterList);

% start execution
fprintf(fid,'%% start the execution');
fprintf(fid,'\r\n');
switch Workflow.Type
    case 'popModel'
        fprintf(fid,sprintf('runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,dataFiles,sensParameterList);'));
    case 'meanModel'
        fprintf(fid,'runMeanModelWorkflow(WSettings,TaskList,MeanModelSet,VPC,dataFiles,sensParameterList,MBS);');

end
fprintf(fid,'\r\n');

fclose(fid);


function    t = getDefinitonTextPopRunset(PopRunSet,iSet)

fn = fieldnames(PopRunSet);

t = sprintf('PopRunSet(%d) = struct(',iSet);
for iFN = 1:length(fn)
    if islogical(PopRunSet(iSet).(fn{iFN})) || isnumeric(PopRunSet(iSet).(fn{iFN}))
        t = sprintf('%s''%s'',%d,',t,fn{iFN},PopRunSet(iSet).(fn{iFN}));
    elseif isstruct(PopRunSet(iSet).(fn{iFN}))
        t = sprintf('%s''%s'',%s,',t,fn{iFN},fn{iFN});
    else        
        t = sprintf('%s''%s'',''%s'',',t,fn{iFN},PopRunSet(iSet).(fn{iFN}));
    end
end
t = sprintf('%s);',t(1:end-1));


return

function [t,cellText] = getDefinitonTextOutputList(OutputList,iO)

fn = fieldnames(OutputList);

t = sprintf('OutputList(%d) = struct(',iO);
cellText = {};
for iFN = 1:length(fn)
    if islogical(OutputList(iO).(fn{iFN})) || isnumeric(OutputList(iO).(fn{iFN}))
        t = sprintf('%s''%s'',%d,',t,fn{iFN},OutputList(iO).(fn{iFN}));
    elseif iscell(OutputList(iO).(fn{iFN}))
        tmp = OutputList(iO).(fn{iFN});
        if isempty(tmp)
            t = sprintf('%s''%s'',[],',t,fn{iFN});
        else
            
            c = sprintf('OutputList(%d).%s = {',iO,fn{iFN});
            for iRow = 1:size(tmp,1)
                c = sprintf('%s''%s'';',c,strjoin(tmp(iRow,:),''','''));
            end
            c = sprintf('%s};',c(1:end-1));
            c = strrep(c,'%','%%');
            cellText{end+1} = c; %#ok<AGROW>
        
            t = sprintf('%s''%s'',[],',t,fn{iFN});
        end
    else        
        t = sprintf('%s''%s'',''%s'',',t,fn{iFN},OutputList(iO).(fn{iFN}));
    end
end
t = sprintf('%s);',t(1:end-1));

t = strrep(t,'%','%%');

return

function  writePopulations(fid,PopRunSet)

% write Definitions of sets of Populations
fprintf(fid,'%% Definitions of sets of Populations');
fprintf(fid,'\r\n');
fprintf(fid,'clear PopRunSet');
fprintf(fid,'\r\n');
for iSet = 1:length(PopRunSet)
    
    for iO = 1:length(PopRunSet(iSet).OutputList)
        [outputListText,cellText] = getDefinitonTextOutputList(PopRunSet(iSet).OutputList,iO);
        fprintf(fid,sprintf('%s',outputListText));
        fprintf(fid,'\r\n');
        for iC = 1:length(cellText)
            fprintf(fid,sprintf('%s',cellText{iC}));
            fprintf(fid,'\r\n');
        end
    end
    
    popRunSetText = getDefinitonTextPopRunset(PopRunSet,iSet);
    fprintf(fid,sprintf('%s',popRunSetText));
    fprintf(fid,'\r\n');
    fprintf(fid,'clear OutputList');
    fprintf(fid,'\r\n');
end
fprintf(fid,'\r\n');

return

function writeMeanModel(fid,MeanModelSet)

% write Definitions of sets of Simulations
fprintf(fid,'%% Definitions of sets of Mean Model simulations');
fprintf(fid,'\r\n'); 
fprintf(fid,'clear MeanModelSet');
fprintf(fid,'\r\n');
for iSet = 1:length(MeanModelSet)
    
    if ~isempty(MeanModelSet(iSet).OutputList)
        for iO = 1:length(MeanModelSet(iSet).OutputList)
            [outputListText,cellText] = getDefinitonTextOutputList(MeanModelSet(iSet).OutputList,iO);
            fprintf(fid,sprintf('%s',outputListText));
            fprintf(fid,'\r\n');
            for iC = 1:length(cellText)
                fprintf(fid,sprintf('%s',cellText{iC}));
                fprintf(fid,'\r\n');
            end
        end
    else
        fprintf(fid,'OutputList = [];');
        fprintf(fid,'\r\n');
    end
    
    popRunSetText = getDefinitonTextMeanModel(MeanModelSet,iSet);
    fprintf(fid,sprintf('%s',popRunSetText));
    fprintf(fid,'\r\n');
    fprintf(fid,'clear OutputList');
    fprintf(fid,'\r\n');
end
fprintf(fid,'\r\n');


function    t = getDefinitonTextMeanModel(MeanModelSet,iSet)

fn = fieldnames(MeanModelSet);

t = sprintf('MeanModelSet(%d) = struct(',iSet);
for iFN = 1:length(fn)
    if islogical(MeanModelSet(iSet).(fn{iFN})) || isnumeric(MeanModelSet(iSet).(fn{iFN}))
        t = sprintf('%s''%s'',%d,',t,fn{iFN},MeanModelSet(iSet).(fn{iFN}));
    elseif isstruct(MeanModelSet(iSet).(fn{iFN}))
        t = sprintf('%s''%s'',%s,',t,fn{iFN},fn{iFN});
    else        
        t = sprintf('%s''%s'',''%s'',',t,fn{iFN},MeanModelSet(iSet).(fn{iFN}));
    end
end
t = sprintf('%s);',t(1:end-1));


return

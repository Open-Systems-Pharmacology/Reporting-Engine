function preparePopulationWorkflow(workflowInputxls,sheet)
% PREPAREPOPULATIONWORKFLOW creates out of xls table workflow script
%
%  Inputs 
%   workflowInputxls (string) name of xls file
%   sheet (string) name of sheet, if empty first sheet is taken

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


if ~exist('sheet','var')
    sheet = 1;
end

[PopRunSet,TaskList,workflowMode,dataFiles,sensParameterList] = readWorkflowInput(workflowInputxls,sheet);


%% add Outputs
for iSet = 1:length(PopRunSet)
    PopRunSet(iSet).OutputList = readOutputXls(PopRunSet(iSet).outputxls,PopRunSet(iSet).outputsheet);
end
PopRunSet = rmfield(PopRunSet,'outputxls');
PopRunSet = rmfield(PopRunSet,'outputsheet');



%%write new worflow.m
fid = fopen('workflow.m','w');
% write header
% write header
writeWorkflowHeader(fid,'popModel',workflowMode)

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

% write Definitions of sets of TaskList
writeTaskList(fid,TaskList);


% write Datafiles
writeDataFiles(fid,dataFiles);

%% write function call VPC
fprintf(fid,'%% get Definition of default plots');
fprintf(fid,'\r\n');
fprintf(fid,'if TaskList.doVPC');
fprintf(fid,'\r\n');
fprintf(fid,'    VPC = getDefaultVPCSettings(WSettings,PopRunSet);');
fprintf(fid,'\r\n');
fprintf(fid,'else');
fprintf(fid,'\r\n');
fprintf(fid,'    VPC = [];');
fprintf(fid,'\r\n');
fprintf(fid,'end');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%% write function call Sensitivity
% write sensitivity Aalysis
writeSensitivityParameterList(fid,sensParameterList);

% start execution
fprintf(fid,'%% start the execution');
fprintf(fid,'\r\n');
fprintf(fid,sprintf('runPopulationWorkflow(WSettings,TaskList,PopRunSet,VPC,dataFiles,sensParameterList);'));
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


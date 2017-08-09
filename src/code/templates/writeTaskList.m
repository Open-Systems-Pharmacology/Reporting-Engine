function writeTaskList(fid,TaskList)
% WRITETASKLIST write Tasklist to workflow.m
%
% writeTaskList(fid,TaskList)
%
% Inputs
%   fid (double) id of file
%   TaskList (structure) list of tasks

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


fprintf(fid,'%% Definitions of TaskList');
fprintf(fid,'\r\n');
fn = fieldnames(TaskList);
t = sprintf('TaskList = struct(');
for iFN = 1:length(fn)
    if islogical(TaskList.(fn{iFN})) || isnumeric(TaskList.(fn{iFN}))
        t = sprintf('%s''%s'',%d,',t,fn{iFN},TaskList.(fn{iFN}));
    else
        t = sprintf('%s''%s'',''%s'',',t,fn{iFN},TaskList.(fn{iFN}));
    end
end
t = sprintf('%s);',t(1:end-1));
fprintf(fid,sprintf('%s',t));
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

return
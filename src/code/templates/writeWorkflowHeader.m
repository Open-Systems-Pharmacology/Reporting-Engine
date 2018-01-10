function writeWorkflowHeader(fid,workflowType,workflowMode,documentation)
% WRITEWORKFLOWHEADER write header of a wrokflow script
%
% writeWorkflowHeader(fid,workflowType,workflowMode)
%
% Inputs
%       fid (double) id of file
%       workflowMode (string) chararcetrises type of workflow
%       workflowMode (string) chararcetrises mode of workflow


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% write header
fprintf(fid,'%% Script to start a workflow');
fprintf(fid,'\r\n');
fprintf(fid,'%% Type: %s',workflowModeToText(workflowType,workflowMode));
fprintf(fid,'\r\n');
for iDoc = 1:size(documentation,1)
    fprintf(fid,'%% %s: %s',documentation{iDoc,1}{1},documentation{iDoc,2}{1});
    fprintf(fid,'\r\n');
end
fprintf(fid,'%% Original author: %s %s',getenv('Username'),datestr(now));
fprintf(fid,'\r\n');
fprintf(fid,'%% \r\n');
fprintf(fid,'%%  HOW TO USE');
fprintf(fid,'\r\n');
fprintf(fid,'%%  this script has to be filed in your working directory together with your input files like the simulation xml');
fprintf(fid,'\r\n');
fprintf(fid,'%%  and the poulation csv');
fprintf(fid,'\r\n');
fprintf(fid,'%%  adjust description for your purpose');
fprintf(fid,'\r\n');
fprintf(fid,'%%  set the matlab directory to your working directory');
fprintf(fid,'\r\n');
fprintf(fid,'%%  start the script');
fprintf(fid,'\r\n');
fprintf(fid,'%% \r\n');
fprintf(fid,'%%  Script is intendended to run with Matlab 2017b ');
fprintf(fid,'\r\n');
fprintf(fid,'%% \r\n');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

% global settings
fprintf(fid,'%% global settings');
fprintf(fid,'\r\n'); 
fprintf(fid,'%% there are globale settings which are used in all functions.');
fprintf(fid,'\r\n');
fprintf(fid,'WSettings = getDefaultWorkflowSettings(''%s'',''%s'');',workflowType,workflowMode);
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');


return

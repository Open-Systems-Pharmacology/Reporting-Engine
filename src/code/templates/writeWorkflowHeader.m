function writeWorkflowHeader(fid)
% WRITEWORKFLOWHEADER write header of a wrokflow script
%
% writeWorkflowHeader(fid)
%
% Inputs
% fid (double) id of file


% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


% write header
fprintf(fid,'%% Script to start a workflow');
fprintf(fid,'\r\n');
fprintf(fid,'%% Purpose:');
fprintf(fid,'\r\n');
fprintf(fid,'%% M&S activity:');
fprintf(fid,'\r\n');
fprintf(fid,'%% Validation level:');
fprintf(fid,'\r\n');
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
fprintf(fid,'%%  ! SPM reporting engine runs with Matlab 2013b (linux cluster)');
fprintf(fid,'\r\n');
fprintf(fid,'%% \r\n');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

% global settings
fprintf(fid,'%% global settings');
fprintf(fid,'\r\n'); 
fprintf(fid,'%% there are globale settings which are used in all functions.');
fprintf(fid,'\r\n');
fprintf(fid,'WSettings = getDefaultWorkflowSettings;');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');


return

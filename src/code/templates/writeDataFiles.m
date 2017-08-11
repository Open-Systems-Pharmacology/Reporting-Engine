function writeDataFiles(fid,dataFiles)
% WRITEDATAFILES write datalist to workflow.m
%
% writeDataFiles(fid,dataFiles)
%
% Inputs
%   fid (double) id of file
%   TaskList (structure) list of tasks

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org


fprintf(fid,'%% List of Nonmemfiles:');
fprintf(fid,'\r\n');
fprintf(fid,'%% set to {}, if no data available');
fprintf(fid,'\r\n');
fprintf(fid,'%%  first column nonmem file, second column dictionary, third column datatype');
fprintf(fid,'\r\n');

if isempty(dataFiles)
    fprintf(fid,'dataFiles = {};');
    fprintf(fid,'\r\n');
else

    for iF = 1:size(dataFiles,1)
        fprintf(fid,'dataFiles = {''%s'',''%s'',''%s''};',dataFiles{iF,1},dataFiles{iF,2},dataFiles{iF,3});
        fprintf(fid,'\r\n');
    end
end
fprintf(fid,'\r\n');

return
function path_out = getGlobalPath(path_in)
%GETGLOBALPATH Support function: Check for existence of the file/folder
% and rectify any mispecification from ../ issues
%
%   getGlobalPath(path_in)
%       path_in (string) path of file or folder    
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%% If path_in is not found, path_out will be same as path in
path_out = fullfile(path_in);

% Check if path is a file
if isfile(path_out)
    PathInfo = dir(path_out);
    path_out = fullfile(PathInfo.folder, PathInfo.name);
% Check if path is a folder
elseif isfolder(path_out)
    PathInfo = dir(path_out);
    path_out = fullfile(PathInfo(1).folder);
% Check for mispecificaion of relative folder
elseif strcmp(path_out(1:3), '..\')
    path_out(1:3) = [];
    % Check if path is a file
    if isfile(path_out)
        PathInfo = dir(path_out);
        path_out = fullfile(PathInfo.folder, PathInfo.name);
    % Check if path is a folder
    elseif isfolder(path_out)
    PathInfo = dir(path_out);
    path_out = fullfile(PathInfo(1).folder);
    else
        path_out = fullfile(path_in);
    end
else
    % If path does not exist
    path_out = fullfile(path_in);
    writeToReportLog('WARNING',[path_in 'not found'],'false');
end
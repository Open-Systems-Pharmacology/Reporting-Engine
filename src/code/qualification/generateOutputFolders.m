function  [SectionsOut, Id] = generateOutputFolders(SectionsIn, REOutputFolder)
%GENERATEOUTPUTFOLDERS Support function: 
% Explore tree of sections and generate folders of sections within REOutputFolder 
% linearize the section structure
%
%   generateOutputFolders(SectionsIn, REOutputFolder)
%       SectionsIn (structure or cell) sections to be created
%       REOutputFolder (string) name of folder where the sections
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%% ---------------------------------------------------% 

% Check for Sections Class, 
% if cell: multiple sections, if struct, only one available
if iscell(SectionsIn)    

% Loop on every section
for i = 1:length(SectionsIn)

    % Get Section Id
    Id(i) = SectionsIn{i}.Id;
    SectionsOut(Id(i)).Id=Id(i);
    SectionsOut(Id(i)).Title=SectionsIn{i}.Title;
    
    % Create folder and add its path in the structure
    path = fullfile(REOutputFolder, SectionsIn{i}.Title);
    SectionsOut(Id(i)).path = path;
    mkdir(path);
    writeToReportLog('INFO',['Section ' SectionsIn{i}.Title ' was successfully created'],'false');
    
    % Copy the content of the section within the folder
    SectionsOut(Id(i)).Content = SectionsIn{i}.Content;
    
    % Check for sub-Sections 
    if isfield(SectionsIn{i}, 'Sections')
    % Repeat the process from this folder
    [SubSectionsOut, SubId] = generateOutputFolders(SectionsIn{i}.Sections, path);
    SectionsOut(SubId)=SubSectionsOut(SubId);
	end
end
% if only one section, then Sections is directly a structure
else
    
    % Get Section Id
    Id = SectionsIn.Id;
    SectionsOut(Id).Id=Id;
    SectionsOut(Id).Title=SectionsIn.Title;
    
    % Create folder and add its path in the structure
    path = fullfile(REOutputFolder, SectionsIn.Title);
    SectionsOut(Id).path = path;
    mkdir(path);
    writeToReportLog('INFO',['Section ' SectionsIn.Title ' was successfully created'],'false');
    
    % Copy the content of the section within the folder
    SectionsOut(Id).Content = SectionsIn.Content;
    
    % Check for sub-Sections 
    if isfield(SectionsIn, 'Sections')
    % Repeat the process from this folder
    [SubSectionsOut, SubId] = generateOutputFolders(SectionsIn.Sections, path);
    SectionsOut(SubId)=SubSectionsOut(SubId);
    end
end

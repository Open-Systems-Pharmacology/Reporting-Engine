function  [SectionsOut, SectionIndex] = generateOutputFolders(SectionsIn, REOutputFolder)
%GENERATEOUTPUTFOLDERS Support function:
% Explore tree of sections and generate folders of sections within REOutputFolder
% linearize the section structure
%
%   generateOutputFolders(SectionsIn, REOutputFolder)
%       SectionsIn (structure or cell) sections to be created
%       REOutputFolder (string) name of folder where the sections
%

% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------------

% Standardize all input sections into cell class
if ~iscell(SectionsIn)
    SectionsIn = {SectionsIn};
end

% Initialize output sections
SectionsOut=[];

% Loop on every section
for i = 1:length(SectionsIn)
    
    % Get Section Id
    SectionIndex=length(SectionsOut)+1;
    SectionsOut(SectionIndex).Id=SectionsIn{i}.Id;
    SectionsOut(SectionIndex).Title=SectionsIn{i}.Title;
    
    % Create folder and add its path in the structure
    Path = fullfile(REOutputFolder, removeForbiddenLetters(SectionsIn{i}.Title));
    SectionsOut(SectionIndex).Path = Path;
    mkdir(Path);
    writeToReportLog('INFO',['Section ' SectionsIn{i}.Title ' was successfully created'],'false');
    
    % Copy the content of the section within the folder
    SectionsOut(SectionIndex).Content = SectionsIn{i}.Content;
    
    % Check for sub-Sections
    if isfield(SectionsIn{i}, 'Sections')
        % Repeat the process from this folder
        [SubSectionsOut, SubSectionIndex] = generateOutputFolders(SectionsIn{i}.Sections, Path);
        for j=1:SubSectionIndex
            SectionsOut(SectionIndex+j)=SubSectionsOut(j);
        end
    end
end

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
    
    try
        mkdir(Path);
        % Create title markdown
        fileID = fopen(fullfile(Path, '_title.md'),'wt');
        fprintf(fileID,'%s\n',SectionsOut(SectionIndex).Title);
        fclose(fileID);
    catch exception
        writeToReportLog('ERROR', [sprintf('Error in Section Id %d: Section Title %s could not be created : \n', SectionsOut(SectionIndex).Id, SectionsOut(SectionIndex).Title), exception.message], 'true', exception);
        rethrow(exception);
    end
    
    % Copy the content of the section within the folder
    SectionsOut(SectionIndex).Content = SectionsIn{i}.Content;
    
    try
        copyfile(SectionsOut(SectionIndex).Content, fullfile(Path, '_content.md'));
    catch exception
        writeToReportLog('ERROR', [sprintf('Error in Section Id %d: %s could not be copied : \n', SectionsOut(SectionIndex).Id, SectionsOut(SectionIndex).Content), exception.message], 'true', exception);
        rethrow(exception);
    end
    
    % Check for sub-Sections
    if isfield(SectionsIn{i}, 'Sections')
        % Repeat the process from this folder
        [SubSectionsOut, SubSectionIndex] = generateOutputFolders(SectionsIn{i}.Sections, Path);
        for j=1:SubSectionIndex
            SectionsOut(SectionIndex+j)=SubSectionsOut(j);
        end
    end
end

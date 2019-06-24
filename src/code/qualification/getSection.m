function  [SectionPath, indexed_items] = getSection(Sections, SectionId)
%GETSECTION Support function:
% Get the path of section of specific section id
%
%   getSection(Sections, SectionId)
%       Sections (array of structures) linearized structure of sections
%       SectionId (integer) ID of section
%
%       SectionPath (string) path of the Section
%       indexed_items (integer) number of files staring with integers
%
% Open Systems Pharmacology Suite;  http://open-systems-pharmacology.org

%---------------------------------------------------

% Fetch in Sections the one with Id sectionID
% If not in Sections ouptut []
SectionPath=[];
indexed_items=0;

for i=1:length(Sections)
    if Sections(i).Id == SectionId
        SectionPath = Sections(i).Path;
    end
end

if ~isempty(SectionPath)
    % Count the number of indexed files in folder (with 3 integers)
    listing=dir(SectionPath);
    for i=1:length(listing)
        if length(listing(i).name)>=3
            if ~isempty(str2num(listing(i).name(1:3)))
                indexed_items=indexed_items+1;
            end
        end
    end
end
function [parPathes,parValues,header] = readPopulationCSV(csvfile)
% READPOPULATIONCSV reads a poulation csv file 
%
% Inputs 
%   - csvfile (string) name of the exported file
% Outputs 
%   - parPathes (cellarray):  pathnames of the exported parameters
%   - parValues (double matrix):  values of the exported parameters and individuals
%   - header (cellarray):  metainfo of the exported file

% Open Systems Pharmacology Suite;  support@systems-biology.com
% Date: 14-July-2017

if ~exist(csvfile,'file')
    error('This file does no exist: %s')
end
    
% read indviduals
fid = fopen(csvfile);
tmp=textscan(fid,'%s','delimiter','\n');
C=tmp{1};
fclose(fid);

% Header
tmp=textscan(C{1},'%s','delimiter',';');
header=tmp{1};
tmp=textscan(C{2},'%s','delimiter',';');
header{2}=tmp{1}{1};

% pathnames
tmp=textscan(C{3},'%s','delimiter',';');
parPathes=tmp{1};


% ParameterValues
parValues=nan(length(C)-3,length(parPathes));
for iC=4:length(C)
    tmp=textscan(C{iC},'%s','delimiter',';');
    v=cellfun(@str2double,tmp{1},'UniformOutput', false);
    parValues(iC-3,:)=cell2mat(v);
end


return
